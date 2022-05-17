library(raster)
library(readxl)
library(data.table)


# HyBIS
tabela <- paste0(gdrive, "/Doutorado/Videos HyBIS/tabela.xlsm")


annotations_raw <- setDT(read_excel(tabela, sheet = "main"))
morphotypes_raw <- setDT(read_excel(tabela, sheet = "morphotypes"))
stations_raw <- setDT(read_excel(tabela, sheet = "stations"))


hybis <- annotations_raw[morphotypes_raw, "code" := .(code), on = .(morphotype)]
hybis <- hybis[stations_raw, date := date, on = .(station)]
hybis[, Time := as.POSIXct(
    paste0(as.character(date, format = "%Y/%m/%d"), as.character(time, format = "%H:%M:%S")),
    format = "%Y/%m/%d %H:%M:%S", tz = "GMT")]

hybis <- hybis[code %in% spps, .(station, Time, Ninds, code)]
colnames(hybis)[1] <- "Dive"


# track
track <- fread(paste0(gdrive, "/Doutorado/Videos HyBIS/track clean/R/track_clean.csv"))

# get environment variables from raster
r.path <- paste0(gdrive, "/Doutorado/images and data/BTM_mosaic/")
r.files <- list.files(r.path, patter = "\\.tif$", full.names = TRUE)
r.hybis <- stack(r.files)

# join track and raster
track[, id := cellFromXY(r.hybis, data.frame(x = Lon, y = Lat))]

hybis <- hybis[track, c("Lon", "Lat", "Depth", "id") := .(Lon, Lat, Depth, id), on = .(Dive, Time)]
rm(annotations_raw, morphotypes_raw, stations_raw)


# call habitat table to remove where is stopped
habitat_raw <- setDT(read_excel(tabela, sheet = "habitat"))
habitat_raw[, c("filename", "comments", "time") := NULL]
colnames(habitat_raw)[1:2] <- c("Dive", "start")
habitat_raw[, end := shift(start, -1)]
habitat_raw <- habitat_raw[!(habitat %in% c("parado", "ignorar"))]
habitat_raw[, id := .I]
habitat_raw <- habitat_raw[, .(Time = seq(start, end, by = 1)), by = .(Dive, habitat, id)]

# subset track when is wathing the floor
track[, Time2 := as.character(Time, format = "%H:%M:%S")]
habitat_raw[, Time2 := as.character(Time, format = "%H:%M:%S")]
track2 <- track[habitat_raw, .(Dive, Time, id, Lon, Lat, Depth), on = .(Dive, Time2)]

track2 <- track2[, c(.(Time = first(Time)), lapply(.SD, median)), by = .(Dive, id), .SDcols = c("Lon", "Lat", "Depth")]

# join the track and annotations
hybis <- dcast(hybis, Dive + id ~ code, value.var = "Ninds", fun.aggregate = sum)
hybis <- hybis[track2, on = .(Dive, id)]
setnafill(hybis, fill = 0, cols = spps)
rm(habitat_raw, track2, tabela)

# extract enviorment variables
envs <- extract(r.hybis, hybis$id)

hybis[, (spps) := lapply(.SD, function(x) factor(x > 0, c(TRUE, FALSE), c("presence", "ausence"))), .SDcols = spps]
hybis <- cbind(hybis, envs)
setorder(hybis, Dive, Time)




# Shinkai

#same order as spps
spps2 <- c("sarostegia_oculata")

# read annotations from EIVA
shinkai <- fread(paste0(gdrive, "/Doutorado/Shinkai/coral_amount_per_image_path.csv"))
colnames(shinkai)[1] <- "image"
shinkai[, c("camera", "Video", "Frame") := tstrsplit(image, "/", keep = 6:8)]


# read outputs from extracting the frames
output1 <- fread(paste0(gdrive, "/EIVA/logs_extracao_frames_Shinkai/CAM1/output.csv"))
output1$camera <- "CAM1"
output2 <- fread(paste0(gdrive, "/EIVA/logs_extracao_frames_Shinkai/CAM2/output.csv"))
output2$camera <- "CAM2"
output <- rbind(output1, output2)
output[, Video := paste0("Dive", substr(Station, 3, 6), "_", sapply(strsplit(output$Video, "_"), `[`, 2))]


# read track data
shinkai <- output[shinkai, on = .(camera, Video, Frame)]
colnames(shinkai)[1] <- "Dive"
shinkai <- shinkai[Dive == "6K1338"]

track <- fread(paste0(gdrive, "/Doutorado/Shinkai/R_track_clean.csv"))
track[, c("Date", "Time") := tstrsplit(Time, " ", keep = 1:2)]
track[, Dive := sub("_", "", Dive)]

shinkai <- track[shinkai, on = .(Dive, Time)]
shinkai[, Time := as.POSIXct(paste(Date, Time), format = "%Y-%m-%d %H:%M:%S", tz = "GMT")]

# filter columns
setnames(shinkai, spps2, spps)
shinkai <- shinkai[camera == "CAM1"] # filtrar apenas a camera 1
shinkai <- shinkai[, .SD, .SDcols = c("Dive", "Time", "Lon", "Lat", "Depth", spps)]


rm(output1, output2, output, track)



# get environment variables from raster
r.path <- paste0(gdrive, "/Doutorado/images and data/BTM_Shinkai/6K_1338/")
r.files <- list.files(r.path, patter = "\\.tif$", full.names = TRUE)
r.shinkai <- stack(r.files)

shinkai[, id := cellFromXY(r.shinkai, data.frame(x = Lon, y = Lat))]


# extrair conteudo das celulas
shinkai <- shinkai[, c(.(Time = first(Time)), lapply(.(Lon = Lon, Lat = Lat, Depth = Depth), median), lapply(.SD, sum)),
                    by = .(Dive, id), .SDcols = spps]
envs <- extract(r.shinkai, shinkai[, id])

shinkai[, (spps) := lapply(.SD, function(x) factor(x > 0, c(TRUE, FALSE), c("presence", "ausence"))), .SDcols = spps]
shinkai <- cbind(shinkai, envs)
setcolorder(shinkai, colnames(hybis))
setorder(shinkai, Dive, Time)

# CHECK
# for NAs in data.frames
if (any(is.na(hybis))) warning("NAs detected in hybis data.frame")

if (any(is.na(shinkai))) warning("NAs detected in shinkai data.frame")

# clear space
rm(envs, r.files, r.path, spps2)
