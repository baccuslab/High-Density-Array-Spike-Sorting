This is an easy example for the bayes optimal online spike sorter (BSS).
The function "E01_simulateToyData" will create an artificial spike template
and a short piece of toy data containing two instances of that template
(spikes). Then, the data is sorted using the BSS.
This example illustrates
a) the way how to simulate a short piece of toy data and the basic usage
   of the BSS sorter.
b) the impact of overlapping spikes of the same template. This is
   important to understand the "dead time" of the spike detection. If a
   neuron is e.g. bursting, single spikes of the same neuron can be very
   close to each other. Can they still be detected?
c) the dependence of the template length (Tf) on the signal to noise
   ratio
d) the dependence of multiple channels on the signal to noise ratio
Remark: c) and d) hold only so nicely if the template is the same on 
        all channels and the noise is white and Gaussian (also accross
        cahnnels).
