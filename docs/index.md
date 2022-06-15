
# 06/15/2022 Update

This week I created a new idealized estuary grid and ran the model for 40 days.
The new grid is straight with a parabolic cross section:

<img src="alpeGrid.png" alt="alpeGrid" width="300"/>
<img src="alpe3D.png" alt="alpe3D" width="300"/>

I ran the model between 01/01/2020 through 02/10/2020 with the following forcing:

|Forcing | Model Input|
|---|---|
|River Flowrate|1000*|
|Tidal Forcing| Includes spring-neap cycle|
|Ocean Forcing| Estuary half fresh at t=0|
|Atmospheric Forcing|Not implemented in this model|

  **Note: Need clarification on units of river flowrate (from river forcing file):*

  ~~~
   Q_mat[:,ii] = 1000 * np.ones(NT) * gri_df.loc[rn, 'isign']
  ~~~

## Results

A salinity cross section plot of the last 10 days is shown below. I noticed that at the end of these 10 days, the estuary becomes more stratified.

<video src="https://user-images.githubusercontent.com/15829099/173406059-3dc31852-9f74-460e-9f1d-61ec65cab31a.mp4" controls="controls" style="max-width: 700px;">
</video>