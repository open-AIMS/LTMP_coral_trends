source('CoralTrends_functions.R')
CoralTrends_checkPackages()

whagbr.n = ML_gClip(whagbr, cbind(c(142,-15.4),c(155,0)))
#whagbr.c = ML_gClip(whagbr, cbind(c(142,-15.4),c(155,-20)))
#whagbr.c = ML_gClip(whagbr, cbind(c(142,-21),c(155,-15.4)))
bb=rbind(c(142,-20.7),
         c(148.7,-20.7),
         c(152,-19.6),
         c(152,-15.4),
         c(142,-15.4))
b.poly=SpatialPolygons(list(Polygons(list(Polygon(bb)),ID=1)))
whagbr.c = ML_gClip(whagbr, b.poly)

#whagbr.s = ML_gClip(whagbr, cbind(c(142,-25),c(155,-20)))
bb=rbind(c(142,-20.7),
         c(148.7,-20.7),
         c(152,-19.6),
         c(155,-19.6),
         c(155,-25),
         c(142,-25))
b.poly=SpatialPolygons(list(Polygons(list(Polygon(bb)),ID=1)))
whagbr.s = ML_gClip(whagbr, b.poly)


data(qld)

save(whagbr, file=paste0(DATA_PATH,'spatial/whagbr.RData'))
save(whagbr.n, file=paste0(DATA_PATH,'spatial/whagbr.n.RData'))
save(whagbr.c, file=paste0(DATA_PATH,'spatial/whagbr.c.RData'))
save(whagbr.s, file=paste0(DATA_PATH,'spatial/whagbr.s.RData'))
save(qld, file=paste0(DATA_PATH,'spatial/qld.RData'))

##Consolidate into one spatial_3Zone
whagbr.n=spChFIDs(whagbr.n,'N')
whagbr.c=spChFIDs(whagbr.c,'C')
whagbr.s=spChFIDs(whagbr.s,'S')
spatial_3Zone = rbind(whagbr.n, whagbr.c, whagbr.s)
#spatial_3Zone = SpatialPointsDataFrame(spatial_3Zone
save(spatial_3Zone, file=paste0(DATA_PATH, 'spatial/spatial_3Zone.RData'))
