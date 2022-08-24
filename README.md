# nCorrect

A toolkit for the correction and normalization of SWC files from neuron morphology experiments.

## Usage

### Parameters

`allAddress`: user-specified text file containing the location of all .swc files (one file path per line)

`imgAddress`: user-specified file location of .tif input image of interest

`outAddress`: user-specified output folder

`a` : Strength of the intensity background term weight (suggested value is `1.00`)

`b` : Strength of the color angle background term weight (suggested value is `0.85`)

`c` : Strength of radius expansion term weight (suggested value is `3.75`)

`T`<sub>`col`</sub> : Threshold of cosine similarity between two colors (suggested value is `0.90`)

`T`<sub>`bkg`</sub> : Threshold for node intensity background (suggested value is `0.47`)

`T`<sub>`min`</sub> : Minimum threshold for lower intensity condition (suggested value is `0.05`)

`T`<sub>`max`</sub> : Maximum threshold for lower intensity condition (suggested value is `0.30`)

`s`<sub>`min`</sub> : Minimum color sum for lower intensity condition (example value is `10000`)

`s`<sub>`max`</sub> : Maximum color sum for lower intensity condition (example value is `85000`)

`r`<sub>`max`</sub> : Maximum node radius (example value is `12`)

### Running

Compile by 

`javac Optimize4D.java`

Run by 

`java Optimize4D allAddress imgAddress outAddress a b c Tcol Tbkg Tmin Tmax smin smax rmax`
