import numpy as np

#Note! One option to resolve the learning rate issue would be to set them all up here. Basically set it so that each parameter moves relatively easily based on the initial values used.


def pinit(batch_size, num_params, num_cells):
    rng = np.random
    p = np.zeros((num_params,num_cells,batch_size))
    lrs = np.zeros((num_params,1,1))
    lr_frac = 0.05

    #All tested starting conditions
    #p[0:4,:] = rng.uniform(0.0, 0.08, size=(4, batch_size)).astype(np.float32) #GSYNs  
    #p[0:4,:] = rng.uniform(0.0, 1.0, size=(4, batch_size)).astype(np.float32) #fPs
    #p[0,:] = rng.uniform(0, 10, size=(1, batch_size)).astype(np.float32) #FRs
    #p[0:4,:] = rng.uniform(-68, -62, size=(4, batch_size)).astype(np.float32) #Els
    #p[0:4,:] = rng.uniform(150, 200, size=(4, batch_size)).astype(np.float32) #Rs
    #p[0:4,:] = rng.uniform(-85, -80, size=(4, batch_size)).astype(np.float32) #Ek
    #p[0:2,:] = rng.uniform(0, 0.3, size=(2, batch_size)).astype(np.float32) #g post IC
    #p[0:4,:] = rng.uniform(30, 80, size=(4, batch_size)).astype(np.float32) #tau
    #p[0:4,:] = rng.uniform(50, 70, size=(4, batch_size)).astype(np.float32) #tauP
    #p[0:4,:] = rng.uniform(5, 20, size=(4, batch_size)).astype(np.float32) #tau_ad
    #p[0:4,:] = rng.uniform(0, 0.01, size=(4, batch_size)).astype(np.float32) #g_inc
    
    #current 4 parameter test
    p[0,:,:] = rng.uniform(0, 0.03, size=(1, num_cells, batch_size)).astype(np.float32) #strf gains
    #p[0,:] = np.ones((1,batch_size)).astype(np.float32)*0.015
    lrs[0] = 0.03*lr_frac

    p[1,:,:] = rng.uniform(0, 0.02, size=(1, num_cells, batch_size)).astype(np.float32) #strf Latencies
    lrs[1] = 0.02*lr_frac

    p[2,:,:] = rng.uniform(0, 0.01, size=(1, num_cells, batch_size)).astype(np.float32) #output adaptation
    lrs[2] = 0.01*lr_frac

    p[3:8,:,:] = rng.uniform(0, 0.2, size=(5, num_cells, batch_size)).astype(np.float32) #GSYNs
    lrs[3:8] = 0.2*lr_frac


    p[8,:,:] = rng.uniform(0, 10, size=(1, num_cells, batch_size)).astype(np.float32) #Absolute refractory period in ms
    lrs[8] = 10*lr_frac

    p[9,:,:] = rng.uniform(0.1, 10, size=(1, num_cells, batch_size)).astype(np.float32) #tanh sharpness
    lrs[9] = 10*lr_frac

    p[10,:,:] = rng.uniform(0, 10, size=(1, num_cells, batch_size)).astype(np.float32) #tanh delay
    lrs[10] = 10*lr_frac

    p[11,:,:] = rng.uniform(0, 1, size=(1, num_cells, batch_size)).astype(np.float32) #Recovery probability return value
    lrs[11] = 1*lr_frac

    # p[4,:] = rng.uniform(5, 20, size=(1, batch_size)).astype(np.float32) #tau_ad
    # lrs[4] = 20*lr_frac
    # p[5,:] = rng.uniform(0.0, 0.01, size=(1, batch_size)).astype(np.float32) #g_inc
    # lrs[5] = 0.01*lr_frac
    #p[6:10,:] = rng.uniform(50, 70, size=(4, batch_size)).astype(np.float32) #tauP
    #lrs[6:10] = 70*lr_frac

    
    
    #min_val = 0
    #max_val = 0.01
    

    #lr = abs(max_val)*0.05


    #---- List of parameters that we definately want to learn ----

    #vvvv   Implemented for single channel  vvvvv
    #FR - Background firing rate 
    #gsyn - strength of synaptic connections
    #fP - Degree of depression   ->  Can be fixed according to gaunonomo (the inhibitory paper that kamal shared a while back)
    #tau - voltage time constant
    #tauP - Depression time constant   -> Immediately dropped MSE to 3.6k? --- Has a very strange warm up period
    #tau_ad - Rate of adaptation  -> Very sensitive to dropping to 0. Adaptation basically disapears at the input if you let it.
    #g_inc - strength of adaptation

    #vvvv  To Validate  vvvvv
    
    #vvvv    To derive   vvvvv
    #t_ref - absolute refractoriness per neuron

    #Number of types of parameters: 9 
    #Number of individual parameters: 29

    #List post meeting of the ones that we want to run for a first pass
    #FR
    #gsyn
    #tau_ad
    #g_inc

    #Maybe taup
    #Maybe gpostIC




    return p,lrs
