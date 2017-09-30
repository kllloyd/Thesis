scatter_fill <- function (x, y, z,xlim=c(min(x),max(x)),ylim=c(min(y),max(y)),zlim=c(min(z),max(z)),
                          nlevels = 20, plot.title, plot.axes, 
                          key.title, key.axes, asp = NA, xaxs = "i", 
                          yaxs = "i", las = 1, 
                          axes = TRUE, frame.plot = axes,xlab='',ylab='', ...) 
{
  #--------------------------------------------------------------------------------------------------------------------------------------------#
  # FROM: https://stackoverflow.com/questions/20127282/r-color-scatterplot-points-by-z-value-with-legend
  # Altered by K Lloyd 2016_08_09
  #--------------------------------------------------------------------------------------------------------------------------------------------#
  mar.orig <- (par.orig <- par(c("mar", "las", "mfrow")))$mar
  on.exit(par(par.orig))
  w <- (3 + mar.orig[2L]) * par("csi") * 2.54
  layout(matrix(c(2, 1), ncol = 2L), widths = c(1, lcm(w)))
  par(las = las)
  mar <- mar.orig
  mar[4L] <- mar[2L]
  mar[2L] <- 1
  par(mar = mar)

# choose colors to interpolate
levels <- seq(zlim[1],zlim[2],length.out = nlevels)
# col <- colorRampPalette(c("red","yellow","dark green"))(nlevels)  
col <- topo.colors(nlevels)  
colz <- col[cut(z,nlevels)]  
#   
plot.new()
plot.window(xlim = c(0, 1), ylim = range(levels), xaxs = "i", yaxs = "i")

rect(0, levels[-length(levels)], 1, levels[-1L],col=col,border=col) 
if (missing(key.axes)) {if (axes){axis(4)}}
       else key.axes
   box()
   if (!missing(key.title)) 
     key.title
   mar <- mar.orig
   mar[4L] <- 1
   par(mar = mar)

   # points
   plot(x,y,type = "n",xaxt='n',yaxt='n',xlab=xlab,ylab=ylab,xlim=xlim,ylim=ylim,bty="n")
   points(x,y,col = colz,xaxt='n',yaxt='n',xlab='',ylab='',bty="n",...)

   ## options to make mapping more customizable

        if (missing(plot.axes)) {
          if (axes) {
            title(main = "", xlab = "", ylab = "")
            Axis(x, side = 1)
            Axis(y, side = 2)
          }
        }
        else plot.axes
        if (frame.plot) 
          box()
        if (missing(plot.title)) 
          title(...)
        else plot.title
        invisible()
 }