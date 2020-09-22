# ROH-DICE

ROH-DICE (<b>R</b>uns <b>o</b>f <b>h</b>omozygous <b>di</b>plotype <b>c</b>luster <b>e</b>numerator) searches for clusters of ROH regions given a genotype panel (in compressed VCF format), minimum length in terms of the number of sites, and the number of samples. 


## Usage:
python ROH-DICE.py [-h] [-i INPUT] [-o OUTPUT] [-w MIN_SAMPLES]
                   [-l MIN_LENGTH] [-L]

<pre>
optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Input file in compressed (.gz) VCF format
  -o OUTPUT, --output OUTPUT
                        Your destination output file.
  -w MIN_SAMPLES, --min_samples MIN_SAMPLES
                        Minimum number of samples in each ROH cluster.
  -l MIN_LENGTH, --min_length MIN_LENGTH
                        Minimum number of sites in each ROH cluster.
  -L, --max_length      Maximize the number of sites in each cluster (by
                        default samples are maximized).


</pre>
