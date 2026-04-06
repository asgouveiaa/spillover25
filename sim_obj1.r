# =============================================================================
# MODELO DE SPILLOVER ZOONÓTICO — NÚCLEO
# Patógenos: E. vogeli (transmissão direta) e Mayaro (transmissão vetorial)
# Estrutura: SIS para animais · Caçadores (H) e Consumidores (C)
# Sazonalidade: precipitação modula K(t), λ(t) e h(t)
# =============================================================================

library(tidyverse)
library(deSolve)
library(ggpubr)

# =============================================================================
# 1. GERAÇÃO DE DADOS DE PRECIPITAÇÃO
# =============================================================================

generate_precipitation <- function(n_years, seed = 42) {
  set.seed(seed)
  total_days <- n_years * 365
  t <- seq(1, total_days)
  
  # Padrão sazonal senoidal + ruído
  base <- 0.5 + 0.4 * sin(2 * pi * (t - 60) / 365)
  noise <- rnorm(total_days, 0, 0.05)
  precip <- pmax(0, pmin(1, base + noise))
  
  data.frame(
    day = t,
    precipitation = base+noise,
    p_normalizado = precip
  )
}

# =============================================================================
# 2. FUNÇÃO DO MODELO — NÚCLEO COMPARTILHADO
# =============================================================================
# Roda E. vogeli OU Mayaro dependendo dos parâmetros:
#   E. vogeli: encounters_day = 0, bite_risk = 0, p_V = 0
#              w_processing > 0, w_eating > 0
#   Mayaro:    encounters_day > 0, bite_risk > 0, p_V > 0
#              w_processing = 0, w_eating = 0

spillover_model_core <- function(t, state, pars, p_data) {
  with(as.list(c(state, pars)), {
    
    # ------------------------------------------------------------------
    # Indexação temporal e precipitação normalizada
    # ------------------------------------------------------------------
    t_int  <- pmax(1, pmin(floor(t), length(p_data)))
    p_norm <- p_data[t_int]
    
    # ------------------------------------------------------------------
    # Sazonalidade climática
    # ------------------------------------------------------------------
    # Capacidade de suporte: reduz na época chuvosa (aglomeração)
    K_effective <- K_max * (K_sazo + (1 - K_sazo) * (1 - p_norm))
    
    # Taxa de transmissão animal: aumenta na época chuvosa
    lambda_effective <- lambda_base * (1 + lambda_sazo * p_norm)
    
    # ------------------------------------------------------------------
    # Dinâmica de caça
    # ------------------------------------------------------------------
    weekday       <- (floor(t) %% 7) + 1
    is_hunt_day   <- ifelse(weekday %in% hunt_days, 1, 0)
    hunting_rate  <- is_hunt_day * h * (1 + h_sazo * p_norm)
    
    # Tempo na floresta por dia de caça (horas/24 → fração do dia)
    time_in_forest <- is_hunt_day * base_time_forest / 24
    
    # ------------------------------------------------------------------
    # Proporção de animais infectados
    # ------------------------------------------------------------------
    total_A    <- A + Ia
    animal_prop <- ifelse(total_A > 0, Ia / total_A, 0)
    
    # ------------------------------------------------------------------
    # DINÂMICA ANIMAL — SIS
    # dA/dt  = crescimento logístico - caça de S - mortalidade - infecção + recuperação
    # dIa/dt = infecção - caça de I  - mortalidade - recuperação
    # ------------------------------------------------------------------
    animal_transmission <- lambda_effective * A * animal_prop
    
    dA  <- r * A * (1 - total_A / K_effective) -
      hunting_rate * A -
      mu * A -
      animal_transmission +
      eta * Ia
    
    dIa <- animal_transmission -
      hunting_rate * Ia -
      mu * Ia -
      eta * Ia
    
    # ------------------------------------------------------------------
    # CAÇADORES — força de infecção
    # Via direta (E. vogeli): processamento e consumo de carne infectada
    # Via vetorial (Mayaro):  picada de vetor infectado na floresta
    # ------------------------------------------------------------------
    direct_risk <- hunting_rate *
      (w_processing * p_P + w_eating * p_E) *
      animal_prop
    
    vector_risk <- time_in_forest *
      encounters_day * bite_risk * p_V *
      animal_prop
    
    hunter_foi <- (direct_risk + vector_risk) * Sh
    
    dSh <- -hunter_foi
    dIh <-  hunter_foi
    
    # ------------------------------------------------------------------
    # CONSUMIDORES — exposição via carne compartilhada
    # ------------------------------------------------------------------
    shared_meat  <- sharing_fraction * hunting_rate * animal_prop
    consumer_foi <- shared_meat *
      (home_processing_risk + home_eating_risk) *
      Sc
    
    dSc <- -consumer_foi
    dIc <-  consumer_foi
    
    # ------------------------------------------------------------------
    # OUTPUT
    # ------------------------------------------------------------------
    list(
      c(dSh = dSh, dIh = dIh,
        dSc = dSc, dIc = dIc,
        dA  = dA,  dIa = dIa),
      
      # Variáveis auxiliares para análise
      hunting_rate     = hunting_rate,
      time_in_forest   = time_in_forest,
      animal_prop      = animal_prop,
      K_effective      = K_effective,
      lambda_effective = lambda_effective,
      direct_risk      = direct_risk,
      vector_risk      = vector_risk,
      prevalence_A     = animal_prop
    )
  })
}

# =============================================================================
# 3. PARÂMETROS — E. VOGELI
# =============================================================================

pars_evogeli <- list(
  # --- Caça ---
  h              = 0.33,       # taxa de caça basal (capturas/caçador/dia de caça)
  h_sazo         = 1.05, #1.2,        # amplificação sazonal da caça
  hunt_days      = c(6, 7),    # dias de caça (6=sábado, 7=domingo)
  base_time_forest = 8,        # horas na floresta por expedição
  
  # --- Dinâmica animal SIS ---
  r              = 0.001,      # taxa de crescimento intrínseco
  K_max          = 50000,      # capacidade de suporte máxima
  K_sazo         = 0.90,       # fração mínima de K na época chuvosa
  mu             = 0.0001,     # mortalidade natural
  eta            = 0.10,       # recuperação animal (SIS)
  lambda_base    = 0.10,       # taxa de transmissão basal entre animais
  lambda_sazo    = 0.30,       # amplificação sazonal da transmissão animal
  
  # --- Spillover direto (E. vogeli) ---
  w_processing   = 1,          # peso da via de processamento
  p_P            = 0.01, #0.05,       # probabilidade de transmissão no processamento
  w_eating       = 1,          # peso da via oral
  p_E            = 0.05, #0.10,       # probabilidade de transmissão oral
  
  # --- Vetorial — zerado para E. vogeli ---
  encounters_day = 0,
  bite_risk      = 0,
  p_V            = 0,
  
  # --- Consumidores ---
  sharing_fraction      = 0.75,   # fração da carne compartilhada
  home_processing_risk  = 0.02,   # risco no processamento doméstico
  home_eating_risk      = 0.10    # risco no consumo doméstico
)

# =============================================================================
# 4. PARÂMETROS — MAYARO
# =============================================================================

pars_mayaro <- list(
  # --- Caça ---
  h              = 0.33,
  h_sazo         = 1.2,
  hunt_days      = c(6, 7),
  base_time_forest = 1,
  
  # --- Dinâmica animal SIS ---
  r              = 0.001,
  K_max          = 50000,
  K_sazo         = 0.90,
  mu             = 0.0001,
  eta            = 0.10,       # recuperação mais lenta para Mayaro
  lambda_base    = 0.12,
  lambda_sazo    = 0.30,
  
  # --- Spillover direto — zerado para Mayaro ---
  w_processing   = 0,
  p_P            = 0,
  w_eating       = 0,
  p_E            = 0,
  
  # --- Vetorial (Mayaro via Haemagogus) ---
  encounters_day = 10,         # picadas por hora de floresta
  bite_risk      = 0.49,       # probabilidade de picada por encontro
  p_V            = 0.15,       # probabilidade de transmissão por picada infectante
  
  # --- Consumidores — sem risco para Mayaro ---
  sharing_fraction      = 0.75,
  home_processing_risk  = 0,
  home_eating_risk      = 0
)

# =============================================================================
# 5. CONDIÇÕES INICIAIS E TEMPO
# =============================================================================

n_years     <- 5
total_days  <- n_years * 365
tempo       <- seq(1, total_days, by = 1)

state_init <- c(
  Sh  = 100,    # caçadores suscetíveis
  Ih  = 0,      # caçadores infectados
  Sc  = 400,    # consumidores suscetíveis
  Ic  = 0,      # consumidores infectados
  A   = 40000,  # animais suscetíveis
  Ia  = 4000    # animais infectados (prevalência inicial ~10%)
)

# =============================================================================
# 6. PRECIPITAÇÃO
# =============================================================================

precipitation_data <- generate_precipitation(n_years)
p_norm_vec <- precipitation_data$p_normalizado

# =============================================================================
# 7. SIMULAÇÕES
# =============================================================================

cat("Rodando modelo E. vogeli...\n")
out_ev <- ode(
  y     = state_init,
  times = tempo,
  func  = spillover_model_core,
  parms = pars_evogeli,
  p_data = p_norm_vec
)

cat("Rodando modelo Mayaro...\n")
out_my <- ode(
  y     = state_init,
  times = tempo,
  func  = spillover_model_core,
  parms = pars_mayaro,
  p_data = p_norm_vec
)

# =============================================================================
# 8. PÓS-PROCESSAMENTO
# =============================================================================

process_output <- function(out, pathogen_name, precip_data) {
  as.data.frame(out) %>%
    left_join(precip_data, by = c("time" = "day")) %>%
    mutate(
      pathogen       = pathogen_name,
      year           = ceiling(time / 365),
      prevalence_H   = ifelse((Sh + Ih) > 0, Ih / (Sh + Ih), 0),
      prevalence_C   = ifelse((Sc + Ic) > 0, Ic / (Sc + Ic), 0),
      prevalence_A   = ifelse((A  + Ia) > 0, Ia / (A  + Ia), 0),
      total_animals  = A + Ia,
      # Garante não-negatividade
      A  = pmax(A,  0),
      Ia = pmax(Ia, 0)
    )
}

results_ev <- process_output(out_ev, "E. vogeli",  precipitation_data)
results_my <- process_output(out_my, "Mayaro",     precipitation_data)
results_all <- bind_rows(results_ev, results_my)

# =============================================================================
# 9. VISUALIZAÇÕES
# =============================================================================

# --- 9.1 Prevalência nos animais ---
p_animal <- results_all %>%
  ggplot(aes(x = time / 365, y = prevalence_A * 100, color = pathogen)) +
  geom_line(linewidth = 0.8) +
  labs(
    title  = "Prevalência no reservatório animal",
    x      = "Tempo (anos)",
    y      = "Prevalência (%)",
    color  = ""
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "top")

# --- 9.2 Proporção de caçadores infectados ---
p_hunters <- results_all %>%
  ggplot(aes(x = time / 365, y = prevalence_H * 100, color = pathogen)) +
  geom_line(linewidth = 0.8) +
  labs(
    title  = "Proporção de caçadores infectados",
    x      = "Tempo (anos)",
    y      = "Proporção infectados (%)",
    color  = ""
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "top")

# --- 9.3 Proporção de consumidores infectados ---
p_consumers <- results_all %>%
  ggplot(aes(x = time / 365, y = prevalence_C * 100, color = pathogen)) +
  geom_line(linewidth = 0.8) +
  labs(
    title  = "Proporção de consumidores infectados",
    x      = "Tempo (anos)",
    y      = "Proporção infectados (%)",
    color  = ""
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "top")

# --- 9.4 Sazonalidade da precipitação ---
p_precip <- precipitation_data %>%
  filter(day <= 365 * n_years) %>%
  ggplot(aes(x = day / 365, y = precipitation)) +
  geom_line(color = "steelblue", linewidth = 0.6) +
  labs(
    title = "Precipitação normalizada (primeiros 3 anos)",
    x     = "Tempo (anos)",
    y     = "Precipitação (normalizada)"
  ) +
  theme_bw(base_size = 12)

# --- 9.5 Painel completo ---
fig_main <- ggarrange(
  p_precip, p_animal,
  p_hunters, p_consumers,
  ncol = 2, nrow = 2,
  common.legend = TRUE,
  legend = "top",
  labels = "auto"
)

print(fig_main)

# =============================================================================
# 10. RESUMO NUMÉRICO
# =============================================================================

summary_results <- results_all %>%
  group_by(pathogen) %>%
  summarise(
    prev_animal_max   = max(prevalence_A, na.rm = TRUE),
    prev_animal_mean  = mean(prevalence_A, na.rm = TRUE),
    prev_hunters_max  = max(prevalence_H, na.rm = TRUE),
    prev_hunters_final = last(prevalence_H),
    prev_consumers_max = max(prevalence_C, na.rm = TRUE),
    animais_final     = last(total_animals),
    .groups = "drop"
  )

print(summary_results)

cat("\nSimulação concluída.\n")
cat("Objetos disponíveis: results_ev, results_my, results_all\n")
cat("Figuras: p_animal, p_hunters, p_consumers, fig_main\n")