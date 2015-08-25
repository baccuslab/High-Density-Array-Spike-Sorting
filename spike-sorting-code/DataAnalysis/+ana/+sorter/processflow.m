1. load enough data to learn spike detector. keep data
2. process buffered data from 1. to detect spikes. and cut them. delete data
3. load and process enough data to cluster spikes and build classifier
4. process all previously detected and cut spikes with classifier. delete spikes
5. process rest.


- loop over files
  - loop over chunks in files
    - sort each chunk
  - assemble one sorting per file and store