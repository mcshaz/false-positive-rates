# false-positive-rates

Explorations of false positives using R

- [False positive rate - Wikipedia](https://en.wikipedia.org/wiki/False_positive_rate)

## Project Status

Currently under casual development.

## Contents

- [advanced-false-positive-risk.R](advanced-false-positive-risk.R)
	- False positive rates, graphing, and phi calculation
- [simple-false-positive-risk.R](simple-false-positive-risk.R)
	- Initial work showing a simple calculation of false positive rates


## What to Expect

### control vs control (p=0.5 vs p=0.5)

typically:
	x-squared < 3.8
	p_value > 0.5


### control vs treatment (p=0.5 vs p=0.8) 

typically:
	x-squared  3.8 to 15
	p_value < 0.5