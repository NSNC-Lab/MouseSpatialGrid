import numpy as np

def adam_update(m,v,p,t,beta1,beta2,lr,eps,out_grad):
   
    #print(p)
    #print(out_grad)
    #print(np.shape(out_grad))

    #Note! Squeezing the channel dimention for now for the single channel work. Couldn't find a way around this. DO NOT SQUEEZE IN FUTURE
    out_grad = np.squeeze(out_grad)

    #print(out_grad)

    #print(np.shape(out_grad))

    #print(np.max(out_grad))
    #print(np.min(out_grad))

    #print('p')
    #print(p)

    #print(lr)

    t += 1
    m = beta1*m + (1-beta1) * out_grad
    v = beta2*v + (1-beta2) * (out_grad**2)

    m_hat = m/(1 - beta1 ** t)
    v_hat = v/(1 - beta2 ** t)

    p = p - lr*m_hat/(np.sqrt(v_hat) + eps)


    #print(lr*m_hat/(np.sqrt(v_hat) + eps))

    #print(out_grad)

    #print(p)

    #print(t)
    #print(m)
    #print(v)

    #print(m_hat)
    #print(v_hat)
    
    #Make sure Gsyns do not go below zero
    #Gsyn & FR & g_inc -- Must be larger than 0
    p[0:5,:] = np.where(p[0:5,:]<0,0,p[0:5,:])
    p[9:13,:] = np.where(p[9:13,:]<0,0,p[9:13,:])
    
    #tau_ad -- kind of sensitive
    p[5:9,:] = np.where(p[5:9,:]<1,1,p[5:9,:])
    p[13:17,:] = np.where(p[13:17,:]<1,1,p[13:17,:])
    #p = np.where(p<0,0,p)

    #p = np.where(p>1,1,p)
    
    #print(p[0][0])

    return m,v,p,t
    