from torch import nn
import torch.nn.functional as F
import torch
from torch.utils.data import Dataset
import numpy as np



class Siamese_cnn(nn.Module):
    def __init__(self, config):
        super(Siamese_cnn, self).__init__()
        self.input_dim = config.input_size
        self.h1_dim = config.h1_dim
        self.h2_dim = config.h2_dim
        self.step=config.cnn_step
        self.k_size = config.cnn_size
        self.pool_size = 2
        self.pool_step = 2
        self.cov = nn.Conv1d(1,1,self.k_size,stride=self.step)
        self.cov_dim = int((self.input_dim-self.k_size)/self.step+1)
        self.pool = nn.MaxPool1d(kernel_size=self.pool_size, stride=self.pool_step)
        self.pool_dim = int((self.cov_dim-(self.pool_size-1)-1)/self.pool_step+1)
        self.encoder = nn.Sequential(
            nn.Linear(self.pool_dim, self.h1_dim),
            nn.Linear(self.h1_dim, self.h2_dim),
        )
    
    def forward(self, s1, s2):
        v1 = self.cov(s1).squeeze(1)
        v1 = self.pool(v1).squeeze(1)
        v2 = self.cov(s2).squeeze(1)
        v2 = self.pool(v2).squeeze(1)
        v1 = self.encoder(v1)
        v2 = self.encoder(v2)
        output = F.cosine_similarity(v1,v2)
        output = 1-torch.clamp(output,1e-6).to(torch.float32)
        
        return output


class myDS(Dataset):
    def __init__(self, df, vec_size,n_thres,mask_na=[]):
        # Assign vocabularies.
        self.input1 = df['s1'].tolist()
        self.input2 = df['s2'].tolist()
        self.na1 = df['na1'].tolist()
        self.na2 = df['na2'].tolist()
        self.r1 = df['r1'].tolist()
        self.r2 = df['r2'].tolist()
        self.name = df['name'].tolist()
        self.size = vec_size
        self.n_thres = n_thres
        self.mask_na = mask_na
        self.r_thres = 15
        self.l_thres = 0.01 
        self.base = 5
    
    def __len__(self):
        return len(self.name)
    
    def __getitem__(self, idx):
        """ initialize """
        na1 = np.array(self.na1[idx].split(','))
        na2 = np.array(self.na2[idx].split(','))
        input1 = np.array([float(i) for i in self.input1[idx].split(',')], dtype='float32')
        input2 = np.array([float(i) for i in self.input2[idx].split(',')], dtype='float32')
        r1 = np.array([int(float(i)) for i in self.r1[idx].split(',')], dtype='float32')
        r2 = np.array([int(float(i)) for i in self.r2[idx].split(',')], dtype='float32')
        name = self.name[idx]
        """ length """
        num1 = len(r1)
        num2 = len(r2)
        """ size and length """
        s1 = input1[:num1]+np.minimum(r1,self.r_thres)
        s2 = input2[:num1]+np.minimum(r2,self.r_thres)
        l1 = input1[num1:]
        l2 = input2[num1:]
        """ length coding """
        l1 = self.length_to_code(l1,thres=self.l_thres,base=self.base)
        l2 = self.length_to_code(l2,thres=self.l_thres,base=self.base)
        """ masking strains """
        mask_index1 = self.mask(na1,num1,self.mask_na)
        mask_index2 = self.mask(na2,num2,self.mask_na)
        s1_ids = s1[mask_index1][:self.size]
        s2_ids = s2[mask_index2][:self.size]
        l1_ids = l1[mask_index1][:self.size]
        l2_ids = l2[mask_index2][:self.size]
        r1_ids = r1[mask_index1][:self.size]
        r2_ids = r2[mask_index2][:self.size]
        """ padding """
        s1_ids = self.pad_array(s1_ids,self.size)
        s2_ids = self.pad_array(s2_ids,self.size)
        l1_ids = self.pad_array(l1_ids,self.size)
        l2_ids = self.pad_array(l2_ids,self.size)
        r1_ids = self.pad_array(r1_ids,self.size)
        r2_ids = self.pad_array(r2_ids,self.size)
        l_effi_flag = (l1_ids  == 0  )+ (l2_ids == 0  )
        l_diff = abs(l1_ids-l2_ids)*l_effi_flag
        l_flag = s1_ids>=s2_ids
        l_flag_inv = np.invert(l_flag)        
        """ Combination """
        n_diff = self.inter_node_diff(r1_ids,r2_ids)
        v1 = s1_ids*n_diff + l_flag*l_diff
        v2 = s2_ids*n_diff + l_flag_inv*l_diff
        """ normalization """
        v1,v2 = self.normalize(v1, v2)
        """ final vector """
        seg1 = np.zeros((1, len(v1)), dtype='float32')
        seg2 = np.zeros((1, len(v2)), dtype='float32')
        seg1[0] = v1
        seg2[0] = v2
        return seg1, seg2, name
    

    def normalize(self,a,b):
        thres = self.n_thres
        max_value = np.maximum(a,b)
        min_value = np.minimum(a,b)
        diff = max_value-min_value
        diff = np.maximum(thres,diff)
        s1 = np.round((a-min_value)/diff,6)
        s2 = np.round((b-min_value)/diff,6)
        s1 = s1 + 0.5
        s2 = s2 + 0.5
        return s1,s2
    
    def mask(self,na,num,mask_na):
        index = [True]*num
        if len(mask_na)==0:
            return index
        else:
            for i in range(num):
                leaf = na[i]
                if leaf in mask_na:
                    index[i]=False
        return index
    
    
    def inter_node_diff(self,r1,r2):
        n_diff = abs(np.minimum(r1,self.r_thres)-np.minimum(r2,self.r_thres))
        return n_diff
    
    
    def pad_array(self,arr, n):
        if len(arr) == 0:
            return np.ones(n, dtype=arr.dtype)*0
        else:
            min_value = np.min(arr[-3:])
            padding_length = n - len(arr)
            if padding_length > 0:
                padding = np.full(padding_length, min_value, dtype=arr.dtype)
                return np.concatenate((arr, padding))
            else:
                return arr

    def length_to_code(self,mat_length,thres=0.01,base=5):
        count = mat_length/thres
        # mat_lcode = np.round(count*base)
        flag = count>1
        mat_lcode = np.round(count*flag*base)
        return mat_lcode
        
    
    
