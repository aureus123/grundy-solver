# An optimization tool that computes the Grundy domination number and the Grundy total domination number

This code contains an implementation of a solver that computes a generalized version of the Grundy domination number, plus other tools to generate instances.
It was partially supported by grants:

  PICT-2016-0410 (ANPCyT)
  PID ING538 (UNR)
  443747/2014-8 (CNPq)
  305264/2016-8 (CNPq)
  PNE 011200061.01.00/16 (FUNCAP/CNPq)

## Files and folders 🔧

grundy.cpp - Source code of the solver
gengraph.cpp - Source code of the random graph generator
genkneser.cpp - Source code of the Kneser graph generator
genweb.pl - Perl script that generates a Web graph
Makefile - self-explained ;)
Set1 - set of random instances (10, 20 and 30 vertices).
Set2 - set of random instances (100 and 200 vertices).
Set3 - set of instances from the DIMACS challenge and its complements (ends in "c"), random
       graphs (25 and 50 vertices) and two graphs representing the city of Buenos Aires.
Set4 - set of Kneser graphs of different sizes (up to 800 vertices).

## Requirements 📋

Use "make" to compile the tools. The solver requires IBM ILOG CPLEX 12.7.

## Examples △

1) To obtain a random graph called "random" of 20 vertices with 50% of edge probability:

        gengraph random.graph 20 50

2) To obtain the Petersen graph:

        genkneser petersen.graph 5 2

3) To obtain the web graph with n = 8 and p = 3 (it'll be called W_8_3.graph):

        genweb.pl 8 3

4) To compute the Grundy domination number of "random" with formulation F4 (and show an
   optimal legal sequence):

        grundy 4 random.graph 1

5) To compute the Grundy total domination number of the Petersen graph:

        grundy 4 petersen.graph 0

## Authors and contact information ✒️

- Manoel Campêlo - UFC, mail: mcampelo@lia.ufc.br
- Daniel Severín - UNR ([**@aureus123**](https://github.com/aureus123)), mail: daniel@fceia.unr.edu.ar

Enjoy! :)
