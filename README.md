# Simulator for "Linear precoding with low-resolution DACs for massive MU-MIMO-OFDM downlink"
(c) 2019 Sven Jacobsson and Christoph Studer;
e-mail: sven.jacobsson@ericsson.com and studer@cornell.edu

### Important information

If you are thinking of contacting us, please do not e-mail the author to ask for download instructions, installation guidelines, or the simulator itself. The simulator itself is well-documented and provides the essential information about the code. Note that we will NOT help to debug user-generated code that was not included in the provided software package. If, however, you notice a bug in our code, please be so kind to contact the Author.

The software package is supplied "as is," without any accompanying support services, maintenance, or future updates. We make no warranties, explicit or implicit, that the software contained in this package is free of error or that it will meet your requirements for any particular application. It should not be relied on for any purpose where incorrect results could result in loss of property, personal injury, liability or whatsoever. If you do use our software for any such purpose, it is at your own risk. The authors disclaim all liability of any kind, either direct or consequential, resulting from your use of these programs.

If you are using the simulator (or parts of it) for a publication, then you *must* cite our paper:

[1] S. Jacobsson, G. Durisi, M. Coldrey, and C. Studer, “Linear precoding with low-resolution DACs for massive MU-MIMO-OFDM downlink,” IEEE Trans. Wireless Commun., to appear.

and clearly mention this in your paper.

### How to start a simulation

Simply run

```sh
precoder_linear_ofdm_sim
```

which starts a simulation for a massive MU-MIMO-OFDM system (with 32 BS antennas, 8 users, 300 occupied subcarriers, and OSR = 1.7) with QPSK modulation using ZF precoding (for the case of 1-bit DACs and infinite resolution). 

The simulator runs with predefined parameters. You can specify your own system and simulation parameters by passing your own "par"-structure (see the simulator for an example). Note that we use default parameters for the considered system configuration; if you want to run the simulation with different parameters, then please refer to the MATLAB code for other parameter settings.

We highly recommend you to execute the code step-by-step (using MATLAB's debug mode) in order to get a detailed understanding of the simulator.

### Version history
- Version 0.1: sven.jacobsson@ericsson.com - simplified/commented code for GitHub
