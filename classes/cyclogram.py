from .utils import mean, peaks

class Cyclogram():
    def __init__(self, heat, time):
        self.heat = heat
        self.time = time
        
        self.Qav, _ = mean(time, heat)
        peaks_index = peaks(self.heat > self.Qav)
        self.peaks_amount = len(peaks_index)
        
        self.Q_peak = []
        self.time_peak = []
        self.time_start = []
        for i in range(self.peaks_amount):
            Q_new, time_new = mean(time, heat, peaks_index[i])
            self.Q_peak.append(Q_new)
            self.time_peak.append(time_new)
            self.time_start.append(time[peaks_index[i][0]])