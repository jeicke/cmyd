function dwelltime=sardwell(range,velocity,lamd,resolution);

range=range*1e3;
dwelltime=lamd*range/(2*velocity*resolution);