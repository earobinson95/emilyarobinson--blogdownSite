---
date: "2021-01-01T00:00:00Z"
external_link: ""
image:
  caption:
  focal_point: Smart
# links:
# - icon: twitter
#   icon_pack: fab
#   name: Follow
#   url: https://twitter.com/georgecushen
slides: 
summary: Method for graphical testing adapted from the NYTimes You Draw It feature.
tags:
- Graphical Testing
- Perception
- D3.js
- r2d3
title: 'You Draw It'
url_code: "https://github.com/srvanderplas/Perception-of-Log-Scales/blob/master/you-draw-it-development/you-draw-it-parameter-selection/main-d3v5.js"
url_pdf: ""
url_slides: ""
url_video: ""
---

One way we can evaluate design choices for data visualizations is through the use of graphical tests. We could ask participants to identify differences in graphs, read information off of a chart accurately, use data to make correct real-world decisions, or
predict the next few observations. All of these types of tests require different levels of use and manipulation of the information presented in the chart. Efforts in the field of graphics have developed graphical testing tools and methods such as the lineup protocol to provide a framework for inferential testing [(Buja et. al,  2009)](https://royalsocietypublishing.org/doi/abs/10.1098/rsta.2009.0120).

In 2015, the New York Times developed a [You Draw it feature](https://www.nytimes.com/interactive/2017/04/14/upshot/drug-overdose-epidemic-you-draw-it.html) where readers are asked to input their own assumptions about various metrics and compare how these assumptions relate to reality. The New York Times team utilizes Data Driven Documents (D3) that allows readers to predict these metrics through the use of drawing a line on their computer screen with their mouse. The goal of this project is to establish ‘You Draw It’, adapted from the New York Times feature, as a tool for graphical testing.

The .gif dislpays an example of a 'You Draw It' task plot. Participants are prompted to "Use your mouse to fill in the trend in the yellow box region." The yellow box region moves along as the participant draws their trend-line until the yellow region disappears. Task plots were created using Data Driven Documents (D3), a JavaScript-based graphing framework that facilitates user interaction. We then integrate this into RShiny using the `r2d3` package.

**Test out the development applet [here](https://emily-robinson.shinyapps.io/you-draw-it-test-app/).**

Future work:

:pencil2: Implement the 'You Draw It' method in non-linear settings.

:chart_with_upwards_trend: Evaluate human ability to extrapolate data from trends.

:cloud: Use the tool to understand beliefs of real data such as climate change trends. 

:package: Develop an R package designed for easy implementation of ‘You Draw It’ task plots.

