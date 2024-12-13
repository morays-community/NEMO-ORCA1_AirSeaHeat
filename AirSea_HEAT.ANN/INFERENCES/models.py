import numpy as np
import tensorflow as tf

model_path = '/lustre/fswork/projects/rech/cli/udp79td/local_libs/morays/NEMO-AirSea_Heat/AirSea_HEAT.ANN/INFERENCES'

#       Utils 
# -----------------
def Is_None(*inputs):
    """ Test presence of at least one None in inputs """
    return any(item is None for item in inputs)

def load_model(path=''):
    """ Load Tensor Flow model from keras-saved neural network. """
    try: # in-repo test or in local deployed config dir
        saved_model = tf.saved_model.load(path+'/TF_ANN')
    except:
        saved_model = tf.saved_model.load('TF_ANN')
    infer = saved_model.signatures["serving_default"]
    return infer

def format_inputs(raw_inputs,tmsk):
    """ Remove land cells and flatten remaining 2D grid. """
    mask_bool = tmsk.flatten() > 0.5
    flat_inputs = raw_inputs.reshape(-1,raw_inputs.shape[2])
    flat_inputs = flat_inputs[mask_bool, :]
    return flat_inputs

def format_outputs(raw_outputs,tmsk):
    """ Build complete grid. Fill land cells with zeros and sea cells with NN results.  """
    final_outputs = np.zeros_like(tmsk.flatten())
    pos = np.where( tmsk.flatten() > 0.5 )[0]
    final_outputs[ pos ] = raw_outputs[:,0]
    return final_outputs.reshape(tmsk.shape[0],tmsk.shape[1],tmsk.shape[2])


#       Main Model Routines
# ------------------------------
# Load model
net = load_model(model_path)

# Predictions    
def airsea_heat_correction(airsea_inputs, tmsk):
    """ Infer air-sea heat flux correction from Storto et al. 2024 Neural Network. """
    if Is_None(airsea_inputs):
        return None
    else:
        global net
        # pack inputs
        airsea_pack = format_inputs(airsea_inputs, tmsk)
        airsea_tensor = tf.convert_to_tensor(airsea_pack.astype('float32'))       

        # predict correction with NN
        preds = net(airsea_tensor)

        # unpack outputs
        correction = preds["dense_4"].numpy()
        return format_outputs(correction,tmsk)
    

if __name__ == '__main__' : 
    test_in = np.random.rand(120,100,24)
    tmsk = np.random.rand(120,100,1)
    test_corr = airsea_heat_correction(test_in, tmsk)
    print(f'Returned correction shape: {test_corr.shape}')
    print(f'Test successful')
