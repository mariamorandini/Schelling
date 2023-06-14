# A statistical physics approach to the Schelling model  of Segregation

Sociological and Economic Model Analysis using Statistical Physics: Schelling and BEG Models

## Table of Contents

- [Introduction](#introduction)
- [Main Results](#main-results)
- [Phase Diagram Analysis](#phase-diagram-analysis)
- [Contents and Simulations]
- [Future Directions](#future-directions)

## Introduction

This repository contains the code and analysis for our study on the analogy between sociological and economic models, specifically the Schelling model, and statistical physics models. Our aim is to gain insights into the phase diagrams of these models and explore their applicability in real-world scenarios.

## Main Results

Through a combination of analytical and numerical analysis, our study demonstrates the validity of the analogy between the sociological and economic model, the Schelling model, and statistical physics models. We have obtained valuable insights by comparing the phase diagrams of these models. 

## 

## Phase Diagram Analysis

To lay the foundation for a detailed analysis of the phase diagrams, we analyzed the Blume Capel model in full generality. This analysis provides a simplified case for understanding the phase behavior. We explored various quantities analogous to physical quantities, such as energy and entropy, to develop a comprehensive physical perspective on the system.

## Contents and Simulations
To validate our phase diagram results and gain further insights, we performed simulations of the models under different parameter combinations. These simulations allowed us to analyze the different types of segregation that occurred within the system. By comparing the simulation results with our analytical findings, we confirmed the accuracy and reliability of our phase diagram analysis.

### class_schelling: 
This class allows to create the instances of the schelling city. Its builtin methods allow to compute in a computationally efficient ways different parameters of the configurations given (i.e. the state of the city). 



### equilibria: 
The main objective of this notebook is to explore the attainment of equilibria under different initial conditions in a Schelling city model. In this context, equilibrium is defined according to Gauvin et al. as follows: "Equilibrium can be understood in two distinct situations: (i) the system no longer evolves (fixed point); (ii) the system reaches a stationary state where the fluctuations of the relevant parameters remain small over a significant number of time steps." In our analysis, we will focus on case (ii) and evaluate these fluctuations by calculating the standard deviation of the model's parameters.





## Future Directions

While our study has provided valuable insights, there are several avenues for further exploration and expansion of the analysis. Some possible directions include:

1. **Inclusion of Chemical Potential:** Consider incorporating a chemical potential to regulate the flow of people into and out of a city. This addition would enhance the realism of the model and provide a more accurate representation of real-world dynamics.

2. **Incorporation of Aging and Economic Possibilities:** Explore the impact of parameters such as aging and economic possibilities on the model. By incorporating these factors, we can better understand the dynamics of societal and economic systems and their interplay with segregation.

3. **Comparison with Ising Model with Non-Uniform Magnetic Field:** Extend our analysis by comparing our model with an Ising model featuring a non-uniform magnetic field. This comparison would help us gain insights into the relationship between price potential in a specific environment and the ease of movement between different locations.

We encourage researchers and enthusiasts to build upon our work and explore these directions for further understanding the dynamics of sociological and economic systems using statistical physics models.


