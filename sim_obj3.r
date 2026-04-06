require(sf)
require(terra)

xmin <- 2
ymin <- 2
xmax <- 5
ymax <- 5

coords <- matrix(c(
  xmin, ymin,
  xmax, ymin,
  xmax, ymax,
  xmin, ymax,
  xmin, ymin
), ncol = 2, byrow = TRUE)

poligono <- sf::st_polygon(list(coords))
poligono_sf <- sf::st_sf(geometry = st_sfc(poligono))

espacamento <- 0.05

grid <- st_make_grid(poligono_sf, cellsize = espacamento, what = "centers")

pontos_sf <- st_sf(geometry = grid)

pontos_dentro <- st_intersection(pontos_sf, poligono_sf)

coords_xy <- st_coordinates(pontos_dentro)
pontos_dentro$x <- coords_xy[,1]
pontos_dentro$y <- coords_xy[,2]

r <- rast(
  xmin = xmin, xmax = xmax,
  ymin = ymin, ymax = ymax,
  resolution = espacamento
)

values(r) <- runif(ncell(r))


r_suave <- focal(r, w = matrix(1, 7, 7), fun = mean)

vals <- extract(r_suave, cbind(pontos_dentro$x, pontos_dentro$y))[,1]

pontos_dentro$suit <- vals

pontos_dentro$suit <- scales::rescale(pontos_dentro$suit, to = c(0,1))

cores <- colorRampPalette(c("blue", "green", "red"))(20)
col_index <- cut(pontos_dentro$suit, breaks = 20, labels = FALSE)

sam <- pontos_sf[sample(nrow(pontos_sf), size = 5), ]
st_crs(sam) <- 4326
sam_utm <- st_transform(sam, crs = 32623)  # ajuste a zona UTM correta
buffer_5km <- st_buffer(sam, dist = 5000)

plot(st_geometry(poligono_sf), col = NA, border = "black")
plot(st_geometry(pontos_dentro), add = TRUE, col = cores[col_index], pch = 15)
plot(st_geometry(buffer_5km), add=TRUE, col = 'grey', main = "Buffer de 5 km")
plot(sam, add = TRUE, col = "black", pch = 16)

# resex

t <- sf::read_sf("./resex_tapajos_arapiuns.kml")

mapbiomas <- terra::rast("./resex2016.tif")
mapbiomas2 <- as.data.frame(mapbiomas,xy=T)

mapbiomas2 <- mapbiomas2 %>% filter(classification_2016 != 0) %>%
  mutate(class = case_when(
    classification_2016 == 3 ~ "Formação Florestal",
    classification_2016 == 6 ~ "Floresta Alagável",
    classification_2016 == 9 ~ "Agropecuária",
    classification_2016 == 11 ~ "Campo Alagado e Área Pantanosa",
    classification_2016 == 21 ~ "Agropecuária",
    classification_2016 == 15 ~ "Agropecuária",
    classification_2016 == 41 ~ "Agropecuária",
    classification_2016 == 29 ~ "Afloramento Rochoso",
    classification_2016 == 39 ~ "Agropecuária",
    classification_2016 == 33 ~ "Rio, Lago e Oceano",
    classification_2016 == 49 ~ "Restinga Arbórea",
    classification_2016 == 5 ~ "Mangue",
    classification_2016 == 40 ~ "Agropecuária",
    classification_2016 == 25 ~ "Outras Áreas não Vegetadas",
    classification_2016 == 24 ~ "Área Urbanizada",
    classification_2016 == 30 ~ "Mineração",
    classification_2016 == 31 ~ "Aquicultura"
  ))

mapbiomas_colors <- c(
  "Formação Florestal" = "#1f8d49", 
  "Floresta Alagável" = "#007785",
  "Mangue" = "#04381d",    
  "Agropecuária" = "#ffefc3",    
  "Agropecuária" = "#ffefc3",   
  "Agropecuária" = "#ffefc3",  
  "Área Urbanizada" = "#d4271e",   
  "Outras Áreas não Vegetadas" = "#db4d4f",   
  "Afloramento Rochoso" = "#ffaa5f",  
  "Mineração" = "#9c0027",    
  "Aquicultura" = "#091077",    
  "Rio, Lago e Oceano" = "#2532e4",    
  "Agropecuária" = "#ffefc3",    
  "Agropecuária" = "#ffefc3",    
  "Agropecuária" = "#ffefc3",    
  "Restinga Arbórea" = "#02d659",
  "Campo Alagado e Área Pantanosa" = "#519799"
)