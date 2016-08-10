import pylab as py
import numpy as np
from spike_sort import extract
import cPickle as pickle

sr = 22050
width = 15
def spike_detect(ar,threshold = .008):
    ''' detect spikes in trace using spike_sort extract function '''
    d = np.array([ar])  # pecularity of spike_sort function requires array of trace arrays ([data_array,])
    dic = {'data' : d,'FS' : sr,'n_contacts' : 1}
    spt = extract.detect_spikes(dic,edge = 'falling', thresh = threshold)
    sp_win = [-1,1]
    if len(spt['data'])>0:  # if spikes detected, then align on peak
        spt = extract.align_spikes(dic,spt, sp_win, type = 'max', resample = 1)
    spike_ar = np.zeros(len(spt['data']))
    for i,spike in enumerate(spt['data']):
        spike_ar[i] = int((spike/1000)*sr)  # convert time from milliseconds to samples
    return(spike_ar)


def get_spike_peaks(raw_data,ar,threshold = .008):
    '''return array of peak heights for ar positions'''
    pos_dict={}
    good_pos=[]
    bad_pos=[]
    for pos in ar:
        pos_dict[pos]={}
        region=pull_spike_region(raw_data,pos,width)
        pos_dict[pos]['max']=region.max()
        pos_dict[pos]['min']=region.min()
        delta=region.max()-region.min()
        if delta < 0.25 * threshold:
            print "shit spike"
            bad_pos.append(pos)
            pos_dict[pos]['type'] = 'shit' 
        else:
            good_pos.append(pos)
            pos_dict[pos]['type']=''
    return(pos_dict,good_pos,bad_pos)

def pull_spike_region(raw_data,pos,width):
    up,down = width, width
    if len(raw_data)-1 < pos+width:
        down = len(raw_data) - 1 - pos
    if pos < width:
        up = pos        
    return(raw_data[pos-up:pos+down])


def remove_duplicates(li):
    new_list = []
    for v in li:
        if v in new_list:
            continue
        else:
            new_list.append(v)
    return(new_list)



def classify_stims_spikes(pos_list,pos_dict,stim_thresh=-.04):
    '''updates pos_dict with type="stim",type="ap" '''
    stim_list=[]
    for pos in pos_list:

        if pos_dict[pos]['min']<stim_thresh:
            pos_dict[pos]['type']='stim'
            stim_list.append(pos)
        else:
            pos_dict[pos]['type']='spike'
    print len(stim_list)
    return(stim_list)


class Spikey_matrix:
    def __init__(self,spike_mat):
        self.low_threshold = 0  # threshold for lowest valid spike
        self.high_threshold = 0 # discard values above 
        self.split_threshold = 0 #
        self.data = np.array(spike_mat)
        
    def onpress(self,event):
        print 'key=%s, x=%d, y=%d, xdata=%f, ydata=%f'%(event.key , event.x, event.y, event.xdata, event.ydata)
        if event.key == 'u':
            self.high_threshold = (int(event.xdata),event.ydata)
            print "high_threshold = %s" % self.high_threshold[0]
        elif event.key == 'b':
            self.low_threshold = (int(event.xdata),event.ydata)
            print "low_threshold = %s" %self.low_threshold[0]
        elif event.key == 't':
            self.split_threshold = (int(event.xdata),event.ydata)
            print "changing threshold"

    def plot_data(self,indicies = 0,ylim = (-.01,.06),xlim = (10,35)):
        if not(indicies):
            indicies = np.arange(len(self.data))
        ''' plot raw, add_values (g), rem_values(r), spikes(b)'''
        fig = py.figure()
        ax = fig.add_subplot(111)
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        [py.plot(self.data[i]) for i in indicies]
        key_cid = fig.canvas.mpl_connect('key_press_event', self.onpress)
        fig.show()
    
    
class Spike_dd_matrix_by_ch:
    '''object for selecting start/end of delay discharge '''
    def __init__(self,raw_data,sort_i = 0):
        self.raw_data = raw_data
        self.num_trials = len(raw_data)
        self.start_pos = [-1] * self.num_trials
        self.end_pos = [-1] * self.num_trials
        self.center_pos = [-1] *self.num_trials
        self.gf_start = [-1] *self.num_trials
        self.gf_end = [-1] *self.num_trials
        self.id_end =  [-1] *self.num_trials
        if isinstance(sort_i,(list,np.ndarray)):
            self.sort_i = sort_i
        else:
            self.sort_i = np.arange(len(raw_data))
        self.ax = None
    def onpress(self,event):
        '''event handling'''
        print 'key=%s, x=%d, y=%d, xdata=%f, ydata=%f axes=%s'%(event.key , event.x, event.y, event.xdata, event.ydata,event.inaxes)
        ax_id = event.inaxes.id
        print ax_id
        if event.key == '[':
            print "ch %s starting at %s"%(ax_id,event.xdata)
            self.start_pos[ax_id] = event.xdata
        elif event.key == ']':
            print "ch %s ending at  %s"%(ax_id,event.xdata)
            self.end_pos[ax_id] = event.xdata
        elif event.key == '.':
            print "ch %s centered at %s" %(ax_id,event.xdata)
            self.center_pos[ax_id] = event.xdata
        elif event.key == '9':
            self.gf_start[ax_id] = event.xdata
            print " GF burst start at %s "%(event.xdata)
        elif event.key == '0':
            print " GF burst end at %s "%(event.xdata)
            self.gf_end[ax_id] = event.xdata
        elif event.key == '6':
            print " ID end at %s"%(event.xdata)
            self.id_end[ax_id] = event.xdata

    def plot_data(self):
        ''' plot raw, add_values (g), rem_values(r), spikes(b)'''
        fig,axes = py.subplots(len(self.raw_data),1,sharex=True,sharey=True,squeeze = False)
        for i,(ax,si) in enumerate(zip(axes,self.sort_i)):
            ax[0].plot(self.raw_data[si])
            ax[0].set_ylim(-0.1,0.1)
            ax[0].id = si
        key_cid = fig.canvas.mpl_connect('key_press_event', self.onpress)
        fig.subplots_adjust(hspace = 0,bottom = 0.05,top = 0.95)
        fig.show()

class Spike_dd_matrix:
    '''object for selecting start/end of delay discharge '''
    def __init__(self,raw_data):
        self.raw_data = raw_data
        self.num_ch = len(raw_data)
        self.start_pos = [-1] * self.num_ch
        self.end_pos = [-1] * self.num_ch
        self.center_pos = [-1] *self.num_ch
        self.gf_start = -1 
        self.gf_end = -1 
        self.ax = None
    def onpress(self,event):
        '''event handling'''
        print 'key=%s, x=%d, y=%d, xdata=%f, ydata=%f axes=%s'%(event.key , event.x, event.y, event.xdata, event.ydata,event.inaxes)
        ax_id = event.inaxes.id
        print ax_id
        if event.key == '[':
            print "ch %s starting at %s"%(ax_id,event.xdata)
            self.start_pos[ax_id] = event.xdata
        elif event.key == ']':
            print "ch %s ending at  %s"%(ax_id,event.xdata)
            self.end_pos[ax_id] = event.xdata
        elif event.key == '.':
            print "ch %s centered at %s" %(ax_id,event.xdata)
            self.center_pos[ax_id] = event.xdata
        elif event.key == '9':
            self.gf_start = event.xdata
            print " GF burst start at %s "%(event.xdata)
        elif event.key == '0':
            print " GF burst end at %s "%(event.xdata)
            self.gf_end = event.xdata


    def plot_data(self):
        ''' plot raw, add_values (g), rem_values(r), spikes(b)'''
        fig,axes = py.subplots(len(self.raw_data),1,sharex=True,squeeze = False)
        for i,ax in enumerate(axes):
            ax.plot(self.raw_data[i])
            ax.set_ylim(-0.1,0.1)
            ax.id = i
        key_cid = fig.canvas.mpl_connect('key_press_event', self.onpress)
        fig.subplots_adjust(hspace = 0,bottom = 0.05,top = 0.95)
        fig.show()
        
    
class Spike_dd:
    '''object for selecting start/end of delay discharge '''
    def __init__(self,raw_data):
        self.raw_data = raw_data
        self.start_pos = None
        self.end_pos = None
        self.center_pos = None

    def onpress(self,event):
        '''event handling'''
        print 'key=%s, x=%d, y=%d, xdata=%f, ydata=%f axes=%s'%(event.key , event.x, event.y, event.xdata, event.ydata,event.inaxes)
        if event.key == '[':
            print "starting at %s"%(event.xdata)
            self.start_pos = event.xdata
        elif event.key == ']':
            print "ending at  %s"%(event.xdata)
            self.end_pos = event.xdata
        elif event.key == 'c':
            print "centered at %s" %(event.xdata)
            self.center_pos = event.xdata
          
    def plot_data(self):
        ''' plot raw, add_values (g), rem_values(r), spikes(b)'''
        fig = py.figure()
        ax = fig.add_subplot(111)
        ax.plot(self.raw_data,'k')
        key_cid = fig.canvas.mpl_connect('key_press_event', self.onpress)
        fig.show()


class Spikey:
    '''object for plotting raw data and obtaining user feedback though mouse/keyboard event handling'''
    def __init__(self,raw_data,threshold,stim_threshold):
        self.raw_data = raw_data
        if threshold:
            self.threshold = threshold
        else:
            self.threshold = 0.01
        if stim_threshold:
            self.stim_threshold = stim_threshold
        else: 
            self.stim_threshold = -0.4
        self.spike_list = []
        self.stim_list = []
        self.a_values = []  # values to add at 
        self.r_values = []
        self.p_values = []  # values for exact placement
        self.bad_list = []
        self.num_neurons = 1



    def onpress(self,event):
        print 'key=%s, x=%d, y=%d, xdata=%f, ydata=%f'%(event.key , event.x, event.y, event.xdata, event.ydata)
        if event.key == 'a':
            self.a_values.append(event.xdata)
        elif event.key == 'd':
            self.r_values.append(event.xdata)
        elif event.key == 't':
            print "changing threshold"
            self.threshold = event.ydata
            #self.a_values.extend(spike_detect(self.raw_data,self.threshold))
        elif event.key == 'b':
            self.stim_threshold = event.ydata
        elif event.key == 'w':
            self.r_values.extend(spike_detect(self.raw_data,event.ydata))
        elif event.key == 'p':
            self.p_values.append(event.ydata)

    def wiggle(self,x,win = 3,max_adjust = 300):
        ''' adjust spikes/stims toward peak '''
        adjusted_values = []
        adjust_list = []
        if x == 'stim':
            adjust_list = self.stim_list
            win = 3
        elif x == 'a':
            adjust_list = self.a_values
            win = 10
            print 'adjusting a_values %s' %len(adjust_list)

        for v in adjust_list:
            new_v = v
            region = self.raw_data[v-win:v+win]
            while(self.raw_data[new_v] < region.max()):
                new_v = new_v - win + region.argmax()
                region = self.raw_data[new_v - win: new_v + win]
            if abs(new_v - v) > max_adjust:
                print 'adjusting too much!!'
                new_v = v
            elif new_v in self.stim_list:
                print 'new value is stim'
                new_v = v
            adjusted_values.append(new_v)
        for i,v in enumerate(adjusted_values):
            adjust_list[i] = v
        


    def preprocess(self):
        '''run spike detection and split spikes and stim artifacts '''
        print "detecting all spikes at %s"%self.threshold
        print "detecting stims at %s"%self.stim_threshold
        all_spikes = spike_detect(self.raw_data,self.threshold)  # detect all spikes at low threshold
        pos_dict,good_list,bad_list=get_spike_peaks(self.raw_data,all_spikes) # sort spikes into good/bad by 
        self.bad_list = bad_list
        self.stim_list = classify_stims_spikes(good_list,pos_dict,self.stim_threshold) 
        self.wiggle('stim')
        remove_stims = self.remove_near_peaks(self.stim_list)
        for stim in remove_stims:
            try:
                self.stim_list.remove(stim)
            except(ValueError):
                print "cannot remove %s" %(stim)
        print '''%s total detected\n%s spikes\n%s rejected spikes\n%s detected stims'''%(len(all_spikes),len(good_list),len(bad_list),len(self.stim_list))
        for pos in pos_dict:
            if pos_dict[pos]['type'] == 'spike':
                self.a_values.append(pos)
       
    
    def remove_near_peaks(self,l,win = 5):
        '''given list of peaks, removes lower of two peaks that are within window of each other'''
        remove_list = []
        ar = np.array(l)
        ar.sort()
        for i,v in enumerate(ar[:-1]):
            if ar[i+1] - v < win:
               if self.raw_data[i+1] > self.raw_data[v]: 
                   remove_list.append(v)
               else:
                   remove_list.append(ar[i+1])
        return(remove_list)




    def add_values(self):
        ''' add values from self.add_values to self.spike_list '''
        for val in self.a_values:
            self.spike_list.append(val)
             #self.spike_list.append(val - 50 + self.raw_data[val - 50: val + 50].argmax())
        self.a_values = []
        self.spike_list.extend(self.p_values)
        
    def remove_values(self):
        ''' remove values in self.remove values from self.spike_list'''
        spike_ar = np.array(self.spike_list)
        to_be_removed = []
        for val in self.r_values:
            i = spike_ar.searchsorted(val)
            closest = 0 
            try:
                if val - spike_ar[i] > spike_ar[i+1] - val:
                    closest = spike_ar[i+1]
                else:
                    closest = spike_ar[i]
                if abs(val - closest) > 100:
                    print 'no value close to %s' % val
                    continue
                else:
                    to_be_removed.append(closest)
            except(IndexError):
                print 'IndexError'
                continue
                
        for rem in to_be_removed:
            print "removing %s"%rem
            try:
                self.spike_list.remove(rem)
            except:
                print "%s not in list" % rem
        self.r_values = []
    def plot_data(self):
        ''' plot raw, add_values (g), rem_values(r), spikes(b)'''
        fig = py.figure()
        ax = fig.add_subplot(111)
        ax.plot(self.raw_data,'k')
        ax.plot(self.stim_list,self.raw_data[self.stim_list],'ko')
        ax.plot(self.spike_list,self.raw_data[self.spike_list],'bo')
        ax.plot(self.a_values,self.raw_data[self.a_values],'go')
        ax.plot(self.r_values,self.raw_data[self.r_values],'ro')
        ax.plot(self.bad_list,self.raw_data[self.bad_list],'rx')
        key_cid = fig.canvas.mpl_connect('key_press_event', self.onpress)
        fig.show()

    def update(self):
        '''set current sorting results to spike_list'''
        print 'updating spike_list with new results'
        self.add_values()
        self.remove_values()
        self.spike_list = remove_duplicates(self.spike_list)
        self.stim_list = remove_duplicates(self.stim_list)

    def clear(self):
        '''remove all spike data'''
        self.a_values = []
        self.r_values = []
        self.p_values = []
        self.spike_list = []

    def subtract(self,sub_ar):
        self.raw_data-=sub_ar




if __name__ == '__main__':
    data = pickle.load(open("sample_data.pkl",'r'))
    x = Spikey(data,.001)

        

    raw_input('done')
