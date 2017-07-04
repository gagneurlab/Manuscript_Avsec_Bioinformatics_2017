from keras.models import Sequential, model_from_config
from keras.layers import Dense, Dropout, Activation, Flatten, Merge
#from keras.layers import Input, merge, LSTM
from keras.layers.normalization import BatchNormalization
from keras.layers.advanced_activations import PReLU
from keras.utils import np_utils, generic_utils
from keras.optimizers import SGD, RMSprop, Adadelta, Adagrad, Adam
from keras.layers import normalization
from keras.layers.convolutional import Convolution2D, MaxPooling2D
from keras.layers import LSTM, Bidirectional
from keras.layers.embeddings import Embedding
from keras.layers.convolutional import Convolution2D, MaxPooling2D, Convolution1D, MaxPooling1D
from keras import regularizers
from keras.callbacks import ModelCheckpoint, EarlyStopping
from keras.constraints import maxnorm

def set_cnn_model(input_dim, input_length):
    nbfilter = 16
    model = Sequential()
    model.add(Convolution1D(input_dim=input_dim, input_length=input_length,
                            nb_filter=nbfilter,
                            filter_length=10,
                            border_mode="valid",
                            # activation="relu",
                            subsample_length=1))
    model.add(Activation('relu'))
    model.add(MaxPooling1D(pool_length=3))

    model.add(Dropout(0.5))

    return model

def get_cnn_network():
    '''
     get_feature = theano.function([origin_model.layers[0].input],origin_model.layers[11].get_output(train=False),allow_input_downcast=False)
    feature = get_feature(data)
    '''
    nbfilter = 16
    print('configure cnn network')

    seq_model = set_cnn_model(4, 107)
    struct_model = set_cnn_model(6, 111)
    # pdb.set_trace()
    model = Sequential()
    model.add(Merge([seq_model, struct_model], mode='concat', concat_axis=1))

    model.add(Bidirectional(LSTM(2 * nbfilter)))

    model.add(Dropout(0.10))

    model.add(Dense(nbfilter * 2, activation='relu'))
    print(model.output_shape)

    return model
