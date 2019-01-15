# Simulator for "Nonlinear precoding for phase-quantized constant-envelope massive MU-MIMO-OFDM"
(c) 2018 Sven Jacobsson and Christoph Studer;
e-mail: sven.jacobsson@ericsson.com and studer@cornell.edu

### Important information

If you are thinking of contacting us, please do not e-mail the author to ask for download instructions, installation guidelines, or the simulator itself. The simulator itself is well-documented and provides the essential information about the code. Note that we will NOT help to debug user-generated code that was not included in the provided software package. If, however, you notice a bug in our code, please be so kind to contact the Author.

The software package is supplied "as is," without any accompanying support services, maintenance, or future updates. We make no warranties, explicit or implicit, that the software contained in this package is free of error or that it will meet your requirements for any particular application. It should not be relied on for any purpose where incorrect results could result in loss of property, personal injury, liability or whatsoever. If you do use our software for any such purpose, it is at your own risk. The authors disclaim all liability of any kind, either direct or consequential, resulting from your use of these programs.

If you are using the simulator (or parts of it) for a publication, then you *must* cite our paper:

[1] S. Jacobsson, O. Castaneda, C. Jeon, G. Durisi, and C. Studer, “Nonlinear precoding for phase-quantized constant-envelope massive MU-MIMO-OFDM,” in IEEE Int. Conf. Telecommunications (ICT), Saint- Malo, France, Jun. 2018, to appear.


and clearly mention this in your paper.

### How to start a simulation

Simply run

```sh
precoder_sim_ofdm
```

which starts a simulation for a massive MU-MIMO-OFDM system (with 128 BS antennas, 16 users, 600 occupied subcarriers, and OSR = 1.7) with 16-QAM modulation using ZF precoding (for 2-phase-bit quantization and infinite resolution) and with the SQUID-OFDM precoder proposed in the paper.

The following precoders are currently supported by the simulator:
  - MRT: maximal-ratio transmission (phase-quantized)
  - MRT_inf: maximal-ratio transmission (infinite resolution)
  - SQUID: squared infinity-norm Douglas-Rachford splitting
  - WF: Wiener filter (phase-quantized)
  - WF_inf: Wiener filter (infinite resolution)
  - ZF: zero-forcing (phase-quantized)
  - ZF_inf: zero-forcing (infinite resolution)

The simulator runs with predefined parameters. You can specify your own system and simulation parameters by passing your own "par"-structure (see the simulator for an example). Note that we use default parameters for the considered system configuration; if you want to run the simulation with different parameters, then please refer to the MATLAB code for other parameter settings.

We highly recommend you to execute the code step-by-step (using MATLAB's debug mode) in order to get a detailed understanding of the simulator.

### Version history
- Version 0.1: sven.jacobsson@ericsson.com - simplified/commented code for GitHub
