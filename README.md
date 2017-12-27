# Matlab code to simulate cell performance of capacitive deionization systems
Capacitive deionization (CDI) is an electrochemical system that can be used to remove salts from water.  

The system works in a similar way as a battery.  When we start to charge the system, ions in the solution can be electrically adsorbed in the electrode, reducing the salt concentration in the water.  

The model accounts for the effect of side reactions, diffuion and convection in porous medium, non-ideal electrosorption theory (which incorporates the effect of co-ion repulsion), and electrode surface charged density.  

Both 1-D and 2-D models have been developed.  

The details of this work can be seen in our [paper](./CDI_pulsation_model_Shang_2017.pdf) or my thesis.

If you use this model in your work, please cite our paper:  

'''
@article{shang2017combined,
  title={A Combined Modeling and Experimental Study Assessing the Impact of Fluid Pulsation on Charge and Energy Efficiency in Capacitive Deionization},
  author={Shang, Xia and Cusick, Roland D and Smith, Kyle C},
  journal={Journal of The Electrochemical Society},
  volume={164},
  number={14},
  pages={E536--E547},
  year={2017},
  publisher={The Electrochemical Society}
}
'''
