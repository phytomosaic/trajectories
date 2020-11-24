######################################################################
#
#    Climate offset strategies (PLOTLY for 3D visualization):
#          emission avoidances and silvicultural interventions
#              required to stabilize community compositions
#
#    Method: trace on NPMR surface in 3+ dimensions
#
#    Rob Smith, phytomosaic@gmail.com, 27 Oct 2020
#
##      GNU General Public License, Version 3.0    ###################

###  Climate offset trajectories:
#      N avoidance and stand silvicultural modification as
#      pathways to climate adaptation
# With CE Ellis and LH Geiser and Jill McMurray
# Use FIA data and the Root airscores model
# 3-D regression surface, ball and arrow trajectories
# Need stand info: hardwood ratio, total BA, canopy cover, etc

### preamble
require(plotly)    # for interactive 3D surface plots
# require(shiny)     # for hosting app

### data preparation
rm(list=ls())
# setwd('/home/rob/prj/trajectories/')
d <- read.csv('./n_west_contour.csv', header=T)
names(d) <- c('Step','x','y','z')
d$z[d$z < 0] <- NA
d <- d[d$x < 7.1,]       # remove a slew of NA at high N dep
`reshape_w` <- function (data, ...) {
        if (!is.data.frame(data))
                stop("must be dataframe")
        if (ncol(data) != 3)
                stop("must have 3-column format")
        x <- data[, 1]
        y <- data[, 2]
        z <- data[, 3]
        fx <- as.factor(x)
        fy <- as.factor(y)
        m <- matrix(NA, nrow = nlevels(fx), ncol = nlevels(fy),
                    dimnames = list(levels(fx), levels(fy)))
        m[cbind(fx, fy)] <- z
        list(x = sort(unique(x)), y = sort(unique(y)), z = m)
}
w <- reshape_w(d[,-1])
x <- w$x
y <- w$y
z <- w$z

# ### EXPERIMENTAL: roughen the surface to see more texture
# set.seed(88)
# z <- z + rnorm(prod(dim(z)), 0, 0.02)

### axis labels setup
# xlb <- expression(Nitrogen~(kg~N~ha^{-1}~y^{-1}))
xlb <- 'N deposition'
ylb <- 'Max Aug Temp (\u00B0C)'
# zlb <- 'Airscores'
zlb <- '<b>Dirty/warm</b> \u2190 Airscores \u2192 <b>Clean/cool</b>'

### function definitions
`make_trace` <- function(vvx, vvy) {
        n  <- 9  # interpolation steps
        vx <- c(sapply(1:(length(vvx)-1), function(i) {seq(vvx[i], vvx[i+1], len=n)}))
        vy <- c(sapply(1:(length(vvy)-1), function(i) {seq(vvy[i], vvy[i+1], len=n)}))
        vz <- sapply(1:length(vx), function(k) { # get target z locations from surface
                z[which.min(abs(x - vx[k])), which.min(abs(y - vy[k]))]})
        isna <- is.na(vz) # smooth right over NA values (TODO...)
        vz[!isna]  <- predict(smooth.spline(vz[!isna] ))$y  # smoothing
        return(list(vx=vx, vy=vy, vz=vz))
}
# # isoline (walk along contour, strategy for NO CHANGE in airscores)
# `make_isoline` <- function(zq=3.7) {
#         n   <- 9
#         dev <- abs(z - zq)
#         iso <- which(dev < quantile(dev, prob=0.01, na.rm=T), arr.ind=T)
#         vvx <- x[iso[,1]]
#         vvy <- y[iso[,2]]
#         vx  <- c(sapply(1:(length(vvx)-1), function(i) {seq(vvx[i], vvx[i+1], len=n)}))
#         vy  <- c(sapply(1:(length(vvy)-1), function(i) {seq(vvy[i], vvy[i+1], len=n)}))
#         isnax <- is.na(vx) # smooth right over NA values (TODO...)
#         isnay <- is.na(vy) # smooth right over NA values (TODO...)
#         vx[!isnax]  <- predict(smooth.spline(vx[!isnax] ))$y  # smoothing
#         vy[!isnay]  <- predict(smooth.spline(vy[!isnay] ))$y  # smoothing
#         data.frame(vx=vx, vy=vy)
# }
# `trans3d` <- function (x, y, z, pmat) {
#         tr <- cbind(x, y, z, 1, deparse.level = 0L) %*% pmat
#         list(x = tr[, 1]/tr[, 4], y = tr[, 2]/tr[, 4])
# }
# `surfcol` <- function (x, ngrid, alpha = 0.9, begin = 0.1, end = 0.85, dir = 1,
#                        pal, ...) {
#         if (missing(pal)) {
#                 pal <- viridis::viridis(n = ngrid, alpha = alpha, begin = begin,
#                                         end = end, direction = dir)
#         }
#         xavg <- (x[-1, -1] + x[-1, -(ncol(x) - 1)] + x[-(nrow(x) -
#                                                                  1), -1] + x[-(nrow(x) - 1), -(ncol(x) - 1)])/4
#         pal[cut(xavg, breaks = length(pal), include.lowest = T)]
# }
# `persp_trace` <- function(vx, vy, vz, theta=-35, ...) {
#         p <- persp(x, y, z, theta = theta, phi = 30,
#                    col = surfcol(z, prod(dim(z)), pal=viridis::viridis(99)),
#                    shade = 0.01, border = '#0000001A',
#                    ticktype = 'detailed', xaxs = 'i',
#                    xlab = xlb, ylab = ylb, zlab = 'Airscores')
#         path <- trans3d(vx, vy, vz, pmat = p)
#         lines(path, col=2, lwd=3)
# }
`plty_trace` <- function(vx, vy, vz, showscale=FALSE, tracecol='red',
                         tracename='tracenm', xlb, ylb, zlb) {
        plot_ly(z=t(z), x=x, y=y, hoverinfo = 'none') %>%
                # add regression surface
                add_surface(showscale=showscale,
                            colorbar = list(
                                    title = 'Airscores'
                                    # tickvals = seq(2,10,by=2),
                                    # ticktext = list('2, Clean/cool','4','6','8',
                                    #                 '10, Dirty/warm')
                            )
                ) %>%
                # add contours to regression surface
                add_surface(
                        contours = list(
                                z = list(show = TRUE, highlight = FALSE,
                                         project=list(z=TRUE), start=0, end=99, size=0.5)),
                        showscale=FALSE) %>%
                # add trace of trajectories to surface
                add_trace(x = vx, y = vy, z = vz+0.05,
                          showlegend=T, name=tracename,
                          type ='scatter3d', mode = 'lines',
                          line = list(color=tracecol, width=4)) %>%
                layout(scene = list(
                                # aspectratio=list(x=0.9, y=1, z=0.8),
                                camera = list(
                                        center = list(x=0.1, y=0.1, z=-0.35),
                                        eye = list(x=-1.75, y=-0.5, z=0)),
                                xaxis = list(title = xlb, type='linear'),
                                yaxis = list(title = ylb, type='linear'),
                                zaxis = list(title = zlb, type='linear')))
}
# ### TEST render
# n <- 7
# tr <- make_trace(seq(2, 3, len=n), seq(22,30,len=n))
# plty_trace(tr$vx, tr$vy, tr$vz, showscale=F,
#            tracecol='red', 'Increase N dep', xlb, ylb, zlb)

# ### plot single SCENARIO
# # target xy locations
# tr <- make_trace(vvx=seq(2, 6, len=7), vvy=c(15, 16.5, 17.5, 20, 30, 32, 32.5))
# # persp
# set_par_mercury(1)
# persp_trace(tr$vx, tr$vy, tr$vz)
# # plotly
# plty_trace(tr$vx, tr$vy, tr$vz, T, '#DF536B', 'tracename', xlb, ylb, zlb)


### set up 5 different SCENARIOS

# # strategy that yields NO CHANGE in airscores (isoline contour walk)
# iso <- make_isoline(zq=3.7)
# iso <- iso[iso$vy > 22 & iso$vy < 30 & iso$vx > 1 & iso$vx < 2,]
# po  <- plty_trace(iso$vx, iso$vy, 3.7, showscale=F,
#                   tracecol='darkblue', 'Iso-trajectory', xlb, ylb, zlb)
# increasing N-dep, 8 degrees warming
n <- 7
tr <- make_trace(seq(2, 3, len=n), seq(22,30,len=n))
p0 <- plty_trace(tr$vx, tr$vy, tr$vz, showscale=FALSE,
                 tracecol='red', 'Increase N dep', xlb, ylb, zlb)
# constant N-dep, 8 degrees warming
tr <- make_trace(rep(2, n), seq(22,30,len=n))
p1 <- plty_trace(tr$vx, tr$vy, tr$vz, showscale=F,
                 tracecol='firebrick', 'Business as usual', xlb, ylb, zlb)
# declining N-dep, 8 degrees warming
tr <- make_trace(seq(2, 1, len=n), seq(22,30,len=n))
p2 <- plty_trace(tr$vx, tr$vy, tr$vz, showscale=F,
                 tracecol='purple', 'Mitigate N dep', xlb, ylb, zlb)
# declining N-dep, 4 degrees warming
tr <- make_trace(seq(2, 1, len=n), seq(22,26,len=n))
p3 <- plty_trace(tr$vx, tr$vy, tr$vz, showscale=F,
                 tracecol='cyan', 'Mitigate N+GHG emissions', xlb, ylb, zlb)
# overplot them
p <- subplot(p0, p1, p2, p3)
# annotate the origin, and make the gradient steeper
a <- list(list(x=tr$vx[1], y=tr$vy[1], z=tr$vz[1],
               arrowsize = 1.5, arrowwidth = 0.95, arrowhead = 3,
               ax = 20, ay=20, textangle = 0, arrowcolor = 'gold',
               text = 'Start',
               font = list(color='gold', size=12),
               xanchor='left', yanchor='bottom'))
p <- p %>%
        layout(scene = list(annotations = a,
                            aspectratio = list(x=0.9, y=1, z=1)
                            )
               ) %>%
        config(displayModeBar = FALSE) %>%
        toWebGL() %>%     # for faster rendering in browser
        partial_bundle()  # for smaller filesize and faster load time


# ### OPTION 1 --- render as interactive plotly object
# p


# ### OPTION 2 --- save it locally as HTML
htmlwidgets::saveWidget(as_widget(p), file = './index.html',
                        title='Mitigation trajectories')
# p$dependencies[[5]]



# ### OPTION 3 --- host as shiny app
# txt_desc <- paste0(
#         'Trajectories for climate and air-quality mitigation. Trajectories ',
#         'describe expected pathways of community compositional change under ',
#         'multiple global-change scenarios. Regression surface is from: '
# )
# txt_note <- paste0(
#         'Patience, takes a moment to load....'
# )
# em_cite <- em(paste0('Geiser, Root, Smith, Jovan, St. Clair, and Dillman (',
#                      'unpublished, 2021).'))
# ui <- fluidPage(
#         verticalLayout(
#                 titlePanel(title = 'Mitigation trajectories',
#                            windowTitle = 'Mitigation trajectories'),
#                 p(txt_desc, em_cite,
#                   'Contact:', a('robert.smith3@usda.gov',
#                                 href='robert.smith3@usda.gov',
#                                 target='_blank')
#                 ),
#                 p(txt_note),
#                 br(),
#                 plotlyOutput('graph', width = '100%', height = '95vh'),
#                 br()
#         )
# )
# server <- function(input, output, session){
#         output$graph <- renderPlotly({ p })
# }
# shinyApp(ui, server)   # <<------





# #############################################################################
# ### zoom in by trimming matrix
# d <- d[d$x < 3 & d$y > 20 & d$y < 32,]
# w <- reshape_w(d[,-1])
# x <- w$x
# y <- w$y
# z <- w$z
# ### strategy that yields NO CHANGE in airscores (isoline contour walk)
# iso <- make_isoline(zq=3.7)
# iso <- iso[iso$vy > 22 & iso$vy < 30 & iso$vx > 1 & iso$vx < 2,]
# po  <- plty_trace(iso$vx, iso$vy, 3.7, showscale=F,
#                   tracecol='darkblue', 'Iso-trajectory', xlb, ylb, zlb)
# # increasing N-dep, 8 degrees warming
# n <- 7
# tr <- make_trace(seq(2, 3, len=n), seq(22,30,len=n))
# p0 <- plty_trace(tr$vx, tr$vy, tr$vz, showscale=T,
#                  tracecol='red', 'N increase', xlb, ylb, zlb)
# # constant N-dep, 8 degrees warming
# tr <- make_trace(rep(2, n), seq(22,30,len=n))
# p1 <- plty_trace(tr$vx, tr$vy, tr$vz, showscale=F,
#                  tracecol='firebrick', 'No change', xlb, ylb, zlb)
# # declining N-dep, 8 degrees warming
# tr <- make_trace(seq(2, 1, len=n), seq(22,30,len=n))
# p2 <- plty_trace(tr$vx, tr$vy, tr$vz, showscale=F,
#                  tracecol='purple', 'N mitigation', xlb, ylb, zlb)
# # declining N-dep, 4 degrees warming
# tr <- make_trace(seq(2, 1, len=n), seq(22,26,len=n))
# p3 <- plty_trace(tr$vx, tr$vy, tr$vz, showscale=F,
#                  tracecol='cyan', 'N+GHG mitigation', xlb, ylb, zlb)
# # overplot them
# p <- subplot(po, p0, p1, p2, p3)
# # annotate the origin, and make the gradient steep
# a <- list(list(x=tr$vx[1], y=tr$vy[1], z=tr$vz[1],
#                arrowsize = 1.5, arrowwidth = 0.95, arrowhead = 3,
#                ax = 10, textangle = 0, arrowcolor = 'gold', text = 'Start',
#                font = list(color='gold', size=12), xanchor='left'))
# p <- p %>% layout(scene = list(annotations = a,
#                                aspectratio = list(x=0.9, y=1, z=1)))
# p
#
# # ### save it locally as HTML
# # getwd()
# # htmlwidgets::saveWidget(as_widget(p), '/home/rob/prj/orise_wildstew/wildstew/www/index.html')
#
#
#
#
#
#
#
#
# ### something like what was in the FARM team prezi:
# # png(file=paste0('./fig/z_08_npmr.png'),
# #     wid=5.0, hei=3.25, unit='in', bg='transparent', res=900)
# require(ggplot2)
# zz <- data.frame(x1 = c(0,0,0.2,0.3,0.4),
#                  y1 = c(22,24,26,28,30),
#                  x2 = c(6.65,6.95,7.05,6.9,6.6),
#                  y2 = c(22,24,26,28,30))
# ggplot(d,aes(x,y,z=z)) +
#      geom_contour_filled(breaks = seq(2,10, by=0.5)) +
#      # scale_colour_viridis(option='D', direction=1, begin=0, end=.9, discrete=T, na.value=NA) +
#      geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2), data=zz, size=1,
#                   inherit.aes = F, color=viridis::inferno(5, begin=0.1, end=0.8)) +
#      xlab(expression(Nitrogen ~ (kg ~ N ~ ha^{-1} ~ y^{-1}))) +
#      ylab('Max Aug Temp (\u00B0C)') +
#      guides(fill=guide_legend(
#           title=expression(Airscores ~ (kg ~ N ~ ha^{-1} ~ y^{-1})), ticks=F,
#           title.position='top', title.hjust=0.5, reverse=T) ) +
#      theme_classic() +
#      theme(plot.background  = element_blank(),
#            panel.background  = element_blank(),
#            legend.title      = element_text(size=10),
#            legend.background = element_blank(),
#            legend.key      = element_blank(),
#            legend.key.size = unit(0.25,'line'),
#            legend.text = element_text(size=6),
#            title = element_text(size=14, face='bold'),
#            # plot.margin=unit(c(2,1,2,0),'mm'),
#            text=element_text(family='Routed Gothic', colour='black'),
#            axis.text=element_text(colour='black')) +
#      scale_x_continuous(breaks=0:7, #limits=c(0,7.5),
#                         expand = expansion(0,0)) +
#      scale_y_continuous(breaks=seq(15,35,by=5), expand = expansion(0,0))
# # dev.off()

####    END    ####
