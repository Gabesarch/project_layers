project_layers

This is a project to develop methods to better identify cortical layers on multi-channel electrode recordings.

+csd is a layer identification toolbox for obtaining superficial, input, and deep layers from multi-channel electrodes.
- get_csd obtains the current source density of the voltage traces. Based on the sink/source reversal point, layer 4 can be identified. 
- get_gamma is a novel method I developed to identify superficial, input, and deep layers from the voltage traces (LFPs) on multi-channel electrodes. Unlike CSD, it is not stimulus dependent. 
