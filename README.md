# factoring
Implementation of some algorithms to factor composite numbers. 
In particular the following algorithms have been implemented:

- Trial Divsion Algorithm
- Fermat's Factoring Algorithm
- Pollard's Rho Factoring Algorithm
- Shanks' Factoring Algorithm


* This repo is in development *


## Synopsis

This project provides the implementation of several factoring algorithms in C++.

## Motivation

This project was started as an assignment for the (DD2440 Advanced Algorithms)[https://www.kth.se/student/kurser/kurs/DD2440?l=en] course at KTH Royal institute of technology. 


## Test

To run the code, first set the value of variable ALGORITHM in line 36 to the desired factoring algorithm. 
Next, compile it writing typing 
`make`

Finally run it and test some sample data.
`./factor < factoring.in`
If you want to test it on specific numbers, you can also run
`python -c "pring NUMBER" | ./factor`
where NUMBER is the number to test.

If you want know the running time, replace `./factor` by `./factor 1`.


## Contributors

If you want to contribute or if you have any question/comment please do not hesitate contacting.


## License

A short snippet describing the license (MIT, Apache, etc.)
