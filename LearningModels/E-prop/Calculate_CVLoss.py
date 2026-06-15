import numpy as np
def calculate(forwards_output, grads,data):
    #Calculate CV of forwards_output and data and then multiply it by the appropriate gradients

    forwards_output = np.squeeze(forwards_output)
    data = np.squeeze(data)

    cv_grads = np.zeros((forwards_output.shape[0],forwards_output.shape[1]))
    cv_loss = 0

    for i in range(forwards_output.shape[0]):
        for j in range(forwards_output.shape[1]):
            data_holder = []
            sim_holder = []
            for k in range(forwards_output.shape[2]):
                data_isis = np.diff(np.flatnonzero(data[i,k,:]))
                sim_isis = np.diff(np.flatnonzero(forwards_output[i,j,k,:]))
                data_holder.append(data_isis)
                sim_holder.append(sim_isis)

            data_holder = [x for x in data_holder if x.size > 0]
            sim_holder = [x for x in sim_holder if x.size > 0]

            if data_holder and sim_holder:
                data_holder = np.concatenate(data_holder)
                sim_holder = np.concatenate(sim_holder)
                cv_grad_val = (sim_holder.std() / sim_holder.mean()) - (data_holder.std() / data_holder.mean())
                cv_loss += cv_grad_val**2
            else:
                cv_grad_val = 0

            if np.isnan(cv_grad_val):
                cv_grad_val = 0
            
            cv_grads[i,j] = cv_grad_val

    print('--------')
    print('cv_loss')
    print(cv_loss)
    print('--------')

    grads[8:12,:,:] = grads[8:12,:,:] * cv_grads[None,:,:,None]

    
    return grads
