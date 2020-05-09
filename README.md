# false-positive-rates

Explorations of false positives using R


## Project Status

Currently under casual development.

## Contents

- [advanced-false-positives.R](advanced-false-positives.R)
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