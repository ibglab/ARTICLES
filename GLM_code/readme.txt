GLM code 

System requirements:
The code runs on Matlab version 2014b.
The optimization toolbox in Matlab is required.
The code does not require any non-standard hardware

Installation:
No installation is needed..

Instructions and Demo:
The main function is "glmMain", which gets 2 inputs: 
1.	spikeTrain - Spike train trials (nTrials X kTimeStamps)
2.	stim - Stimulus of length k
3.	trainPercent – percent of time to use for train data (a number between 0-1). 
Default: 0.8
4.	postSpikeFilterLength – length (in ms) of the post spike filter
Default: 100
5.	stimFilterLength- length (in ms) of the stimulus filter
Default: 60
6.	baseBinWidth – width of each boxcar basis function of the filters
Default: 4

The function generates the GLM, and validates it on the test data.
The outputs of the main function:
1.	stimFilter – the stimulus filter
2.	postSpikeFilter – the post spike filter
3.	rateBias – the bias term
4.	psthCorr- the PCC between the PSTH's of the test data and the simulated spike trains.
For demo, load and run the code on the example data: exampleData.mat (the spike train data) and exampleStim.mat (example stimulus).
The expected output for running the demo with the default values:
1. stimFilter: [0.0076, 0.0073, 0.0067 ,....]
2. postSpikeFilter: [-15.6085  -15.6339  -15.6943, ....]
3. rateBias:    -5.8602
4. psthCorr:     0.9333


