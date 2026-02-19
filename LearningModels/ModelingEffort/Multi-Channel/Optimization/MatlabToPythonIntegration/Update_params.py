import numpy as np

def adam_update(m,v,p,t,beta1,beta2,lr,eps,out_grad):
   
    print(out_grad)

    t += 1
    m = beta1*m + (1-beta1) * out_grad
    v = beta2*v + (1-beta2) * (out_grad**2)

    m_hat = m/(1 - beta1 ** t)
    v_hat = v/(1 - beta2 ** t)

    p = p - lr*m_hat/(np.sqrt(v_hat) + eps)

    return m,v,p,t
    