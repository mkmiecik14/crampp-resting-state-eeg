# EEG Topography Tools
# Matt Kmiecik
# Started 3 June 2022

# HEAVILY inspired by Matt Craddock's code (i.e., credit should go to him):
# https://www.mattcraddock.com/blog/2017/02/25/erp-visualization-creating-topographical-scalp-maps-part-1/
# I also used a similar script in the CRAMPP SSVEP project:
# https://github.com/mkmiecik14/ssvep/blob/main/topo_tools.R

# loads packages and prepares workspace
source("r-prep.R")

# Electrode Locations
# Loads in electrode positions and converts from polar to cartesian
# note: use little x and little y for the coordinates, not X and Y from MATLAB
elec_locs <- 
  read_excel(path = "../data/ss-info.xlsx", sheet = "elec") %>%
  mutate(
    radianTheta = pi/180*theta, 
    x = radius*sin(radianTheta), 
    y = radius*cos(radianTheta)
  ) %>%
  # removes the ' in the labels from matlab to create numeric values for joins
  mutate(labels = gsub("'","",labels))

# theme_topo()
# Plotting tools for the topo map with head shape:
theme_topo <- 
  function(base_size = 12){
    theme_bw(base_size = base_size) %+replace%
      theme(
        rect = element_blank(),
        line = element_blank(), 
        axis.text = element_blank(),
        axis.title = element_blank()
      )
    }

# circleFun()
# function for drawing head for topo plots
circleFun <- 
  function(center = c(0,0), diameter = 1, npoints = 100){
    r = diameter / 2
    tt <- seq(0, 2*pi, length.out = npoints)
    xx <- center[1] + r * cos(tt)
    yy <- center[2] + r * sin(tt)
    return(data.frame(x = xx, y = yy))
    }

# Headshape coordinates are determined
headShape <- 
  circleFun(
    c(0, 0), 
    diameter = round(max(elec_locs$x)), 
    npoints = 100
  )

# Nose coordinates
nose <- data.frame(x = c(-0.075, 0, .075), y = c(.495, .575, .495))

# Produces blank topography map with electrodes
ggplot(headShape, aes(x, y)) +
  geom_path() +
  geom_text(data = elec_locs, aes(x, y, label = labels)) +
  geom_line(data = nose, aes(x, y, z = NULL)) +
  theme_topo() +
  coord_equal()

# Function for interpolating
# data is a data frame containing electrode coordinates at x and y
# currently is accepting dB
# usage topo_terp(data = data, dv = "dB", gridRes = 67)
topo_interp <- 
  function(data = data, dv = dv, gridRes = 67){
    meas <- data[[dv]] # allows for assignment of different dependent variables
    tmp_topo <- 
      with(
        data,
        akima::interp(
          x = x, 
          y = y, 
          z = meas, # let's see if this works
          xo = seq(min(x)*2, max(x)*2, length = gridRes),
          yo = seq(min(y)*2, max(y)*2, length = gridRes),
          linear = FALSE,
          extrap = TRUE
        )
      )
    # Creating a matrix that is x by y filled with z
    interp_topo <- data.frame(x = tmp_topo$x, tmp_topo$z)
    names(interp_topo)[1:length(tmp_topo$y) + 1] <- tmp_topo$y
    interp_topo <- interp_topo %>% gather(key = y, value = !!dv, -x, convert = TRUE)
    
    # mark grid elements that are outside of the plotting circle
    # changing the number after the less than sign will modulate the size of the
    # topo drawn
    interp_topo$incircle <- sqrt(interp_topo$x^2 + interp_topo$y^2) < .85 # original value was .7
    interp_topo <- interp_topo[interp_topo$incircle,] #remove the elements outside the circle
    return(interp_topo)
}

# Function for plotting EEG scalp topographies
topo_plot <- 
  function(
    orig_data = orig_data, # original data
    interp_data = interp_data, # interpolated data
    dv = dv, # name of the column that has the values to plot
    maskRing = circleFun(diameter = 1.7), # creates the mask
    contour_alpha = 1/3, # alpha level of contour lines
    contour_color = "black", # color of contour lines
    headshape_size = .25, # headshape size
    electrode_size = 1, # size of electrode points
    nose_size = .25, # size of nose shape
    bwidth = .5, # width of colorbar
    bheight = .1, # height of colorbar
    color_pal =  brewer.pal(n = 9, "Purples"),
    color_pal_limits = c(-30, -10),
    color_pal_breaks = c(-30, -20, -10),
    legend_name = "DV",
    elec_shape_col = NULL,
    elec_shapes = NULL
    ){
  
    # # looks to see if elec_shape_col has an entry
    # if(is.null(elec_shape_col)){
    #   elec_shape_plot <- geom_point(data = orig_data, aes(x, y), size = electrode_size)
    #   elec_shapes_manual <- scale_shape_manual(values = c(19))
    # }
    # else{
    #   elec_shape_plot <- geom_point(data = orig_data, aes(x, y, shape = {{elec_shape_col}}), size = electrode_size)
    #   elec_shapes_manual <- scale_shape_manual(values = elec_shapes)
    # }
    
    plot <- 
      ggplot(interp_data, aes(x = x, y = y, fill = {{dv}})) +
      coord_equal() + # equalizes coordinates
      geom_raster(interpolate = TRUE) + # basis of the topo
      stat_contour(aes(z = {{dv}}), colour = contour_color, alpha = contour_alpha) + # contour lines
      # creates a mask to remove points outside defined circle
      geom_path(
        data = maskRing,
        aes(x, y, z = NULL, fill = NULL),
        colour = "white",
        size = 6
      ) +
      theme_topo() + # topo theme is added (white background etc.)
      # plots headshape
      geom_path(data = headShape, aes(x, y, z = NULL, fill = NULL), size = headshape_size) +
      # plots nose
      geom_path(data = nose, aes(x, y, z = NULL, fill = NULL), size = nose_size) +
      # colors here
      # note: oob = squish forces everything outside the colour limits to equal
      # nearest colour boundary (i.e., below min colours = min colour)
      scale_fill_gradientn(
        colours = color_pal,
        # may want to hard code these for comparison across groups
        limits = color_pal_limits, # these should be determined from the uninterpolated data
        breaks = color_pal_breaks, labels = color_pal_breaks,
        guide = "colourbar",
        oob = squish,
        name = legend_name
      ) + 
      geom_point(data = orig_data, aes(x, y, shape = {{elec_shape_col}}), size = electrode_size) +
      scale_shape_manual(values = elec_shapes) +
      guides(
        shape = "none", 
        fill = guide_colourbar(
          title.position = "top", 
          title.hjust = 0.5, 
          frame.colour = "black", 
          ticks.colour = "black", 
          barwidth = unit(bwidth, "in"),
          barheight = unit(bheight, "in")
        )
      ) +
      theme(legend.position = "bottom")
    # ANY FACETTING AND/OR SHAPE ADJUSTMENT MUST BE DONE OUTSIDE
    return(plot)
}



