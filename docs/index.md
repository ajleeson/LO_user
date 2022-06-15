# 06/21/2022 Update

### Weekly Todos
 [ ] Create mooring extraction and make a superplot to see tides


# 06/15/2022 Update

This week I created a new idealized estuary grid and ran the model for 40 days.
The new grid is straight with a parabolic cross section:

<img src="alpeGrid.png" alt="alpeGrid" width="300"/>
<img src="alpe3D.png" alt="alpe3D" width="300"/>

I ran the model between 01/01/2020 through 02/10/2020 with the following forcing:

|Forcing | Model Input|
|---|---|
|River Flowrate|1000 cms|
|Tidal Forcing| Includes spring-neap cycle|
|Ocean Forcing| Estuary half fresh at t=0|
|Atmospheric Forcing|Not implemented in this model|

## Results

A salinity cross section plot of the last 10 days is shown below. I noticed that at the end of these 10 days, the estuary becomes more stratified.

<video src="https://user-images.githubusercontent.com/15829099/173406059-3dc31852-9f74-460e-9f1d-61ec65cab31a.mp4" controls="controls" style="max-width: 650px;">
</video>

### Differences between Spring & Neap Tide

Now let's compare stratification during spring and neap tide.

Based on the video above, max tidal amplitudes appear to occur on February 4th. Neap tides likely occur one week later (February 11). Since I only ran the model through February 10, I will look at the estuary on February 10th instead of the 11th. A snapshot the estuary cross section on February 4th and 10th are shown below:

**Spring Tide (Feb 4th)**

<img src="alpeSpring.png" alt="alpeGrid" width="600"/>

**Neap Tide (Feb 10th)**

<img src="alpeNeap.png" alt="alpeNeap" width="600"/>

Based on the plots above, the estuary appears to be more stratified during neap tide compared to spring tide.

### Ongoing Work

Ideally, I would like to make a movie that plots a timeseries of the surface height at a point in the estuary. This way we can visualize the exact dates of spring and neap tides. My attempt at such a movie is shown below (tested for one day, Feb 10th). <span style="color:red">
  I am having trouble with datapoints from previous history files being overwritten by newer history files. What is the best way of making a timeseries from consecutive history files?
</span>

<video src="https://user-images.githubusercontent.com/15829099/173740184-6252fa72-f6a8-4dc2-b053-5898dc5de24b.mp4" controls="controls" style="max-width: 650px;">
</video>