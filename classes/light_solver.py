import numpy as np
import pandas as pd
from scipy import linalg
from scipy import optimize
import itertools

sigma = np.float64(5.67*10**(-8))
earth_radius = np.float64(6371000)
earth_mass = np.float64(5.976 * 10**24)
G = np.float64(6.67 * 10**(-11))
albedo_coef = 0.3
earth_radiation = 240
qs = np.float64(1367)

class Solver:
    def __init__(self, width, Lx, Ly, Lz, c, p, R, orbit_height, cyclograms) -> None:
        self.cyclograms = cyclograms

        self.width = 0.005
        self.Lx = Lx
        self.Ly = Ly
        self.Lz = Lz
        self.R = R * np.ones(shape=(6,6))
        for i in range(6):
            self.R[i][i] = 10**46
        for i in range(3):
            self.R[2*i+1][2*i] = 10**46
            self.R[2*i][2*i+1] = 10**46

        self.areas = np.array(
            [
                Ly*Lz, Ly*Lz,
                Lx*Lz, Lx*Lz,
                Lx*Ly, Lx*Ly
            ]
        )
        self.heat_capacity = c*p*width*self.areas[:]
        
        self.radius = orbit_height*1000 + earth_radius # в метры
        velocity = np.sqrt(G * earth_mass / self.radius) # скорость 
        length = 2*np.pi*self.radius #длина орбиты
        self.period = length/velocity #период полного оборота
        self.angle_sunset = np.arccos(earth_radius/(self.radius))
        self.time_day = self.period * (np.pi/2+self.angle_sunset)/np.pi
        
        self.labels = ["Зенит", "Надир", "Перпенд 1", "Перпенд 2", "Теневая", "Солнечная"]
        
        self.Qav_As = np.zeros(6)
        self.Qav_e = np.zeros(6)
        self.Qmax_As = np.zeros(6)
        self.Qmax_e = np.zeros(6)
        self.Qmin_As = np.zeros(6)  
        self.Qmin_e = np.zeros(6)     
        self.for_12hour_orbit()
        
        self.m = np.zeros(self.R.shape)
        
        self.x0 = 300*np.ones(6)
        self.lambdas = [None]*6
        self.n = 0

    def make_lambda(self, Qav_curr, i, e):
        R_reverse_sum = np.sum(1/self.R[i])
        product_rad = sigma*self.areas[i]*e
        R_reverse = self.R[i]
        result = lambda x: (Qav_curr + np.sum(x/R_reverse)) / (product_rad*x[i]**3+R_reverse_sum) - x[i]
        return result
    
    def func(self, x):
        return list(map(lambda f: f(x), self.lambdas))
        
    def temperature_diapozone(self, As, e, n):
        if self.n % 1000 == 0:
            print(f"Accurate calculation: {self.n} of {n}")
        self.n += 1
        
        result = np.ones(0)
        for i in range(len(self.cyclograms)):
            Qav_iter = self.cyclograms[i].Qav / 6
            Qav = As[:]*self.Qav_As[:] + e[:]*self.Qav_e[:] + Qav_iter
            
            for j in range(6):
                self.lambdas[j] = self.make_lambda(Qav[j], j, e[j])
            
            Tav = optimize.root(self.func, x0=self.x0).x
            
            #расчёт максимума
            Qex = As[:]*self.Qmax_As[:] + e[:]*self.Qmax_e[:] + Qav_iter
            time_day_calc = self.time_day/2
            
            self.m = -1/self.R
            self.m[np.diag_indices_from(self.m)] = np.sum(1/self.R, axis=1) + self.heat_capacity[:]/time_day_calc        
            d = Qex[:] - Qav[:] + Tav[:]*self.heat_capacity[:]/time_day_calc
            Tmax = linalg.solve(self.m, d, overwrite_a=True, overwrite_b=True)
            
            #расчёт минимума
            Qex = As[:]*self.Qmin_As[:] + e[:]*self.Qmin_e[:] + Qav_iter
            time_night_calc = (self.period - self.time_day)/2

            self.m = -1/self.R
            self.m[np.diag_indices_from(self.m)] = np.sum(1/self.R, axis=1) + self.heat_capacity[:]/time_night_calc
            d = Qex[:] - Qav[:] + Tav[:]*self.heat_capacity[:]/time_night_calc
            Tmin = linalg.solve(self.m, d, overwrite_a=True, overwrite_b=True)
            
            iter_result = np.concatenate([Tmin,Tav,Tmax])
            for j in range(self.cyclograms[i].peaks_amount):
                Qex = None
                Q = self.cyclograms[i].Q_peak[j]
                time = self.cyclograms[i].time_peak[j]
                if self.cyclograms[i].time_start[j] < self.time_day:
                    Qex = As[:]*self.Qmax_As[:] + e[:]*self.Qmax_e[:] + (Q / 6)
                else:
                    Qex = As[:]*self.Qmin_As[:] + e[:]*self.Qmin_e[:] + self.cyclograms[i].Q_peak[j]
                
                self.m = -1/self.R
                self.m[np.diag_indices_from(self.m)] = np.sum(1/self.R, axis=1) + self.heat_capacity[:]/time
                d = Qex[:] - Qav[:] + Tav[:]*self.heat_capacity[:]/time_night_calc
                T = linalg.solve(self.m, d, overwrite_a=True, overwrite_b=True)
                iter_result = np.concatenate([iter_result, T])
            
            result = np.concatenate([result, iter_result])
            
        return result
    
    def for_12hour_orbit(self):
        arg = np.sqrt(1-(earth_radius/self.radius)**2)
        angle_coef = 1/(2*np.pi) * (np.pi - 2*np.arcsin(arg) - np.sin(np.arcsin(arg)))
        print(angle_coef)
        #средние
        self.Qav_As[0] = self.areas[0]*qs/np.pi #зенит
        self.Qav_As[1] = self.areas[1]*(1-np.cos(self.angle_sunset))*qs/np.pi #надир после 6 часов
        self.Qav_As[5] = self.areas[5]*(1+np.sin(self.angle_sunset))*qs/np.pi #после до 6 часов
        self.Qav_As[1] += self.areas[1]*qs*albedo_coef/np.pi
        self.Qav_As[2:] = self.Qav_As[2:]+self.areas[2:]*qs*albedo_coef*angle_coef/np.pi
        
        self.Qav_e[1] = self.areas[1]*earth_radiation*earth_radius**2/self.radius**2
        self.Qav_e[2:] = self.Qav_e[2:] + self.areas[2:] * earth_radiation * angle_coef
        
        #максимальные
        self.Qmax_As[0] = qs * self.areas[0] * np.sqrt(2)/2
        self.Qmax_As[5] = qs * self.areas[1] * np.sqrt(2.5)/2
        self.Qmax_As[1] = 2*albedo_coef*qs*self.areas[1]
        self.Qmax_As[2:] = self.Qmax_As[2:]+self.areas[2:]*qs*albedo_coef*angle_coef
        
        self.Qmax_e[1] = self.areas[1]*earth_radiation*earth_radius**2/self.radius**2
        self.Qmax_e[2:] = self.Qmax_e[2:] + self.areas[2:] * earth_radiation * angle_coef
        #минимальные
        self.Qmin_e[1] = self.areas[1]*earth_radiation*earth_radius**2/self.radius**2
        self.Qmin_e[2:] = self.Qmin_e[2:] + self.areas[2:] * earth_radiation * angle_coef
    
    def construct_df(self, As_list, e_list, Tmin, Tmax):
        result_As = []
        result_e = []
        counter = 0
        n = (len(As_list) * len(e_list))**6
        
        for comb in itertools.product(As_list, e_list, repeat=6):
            if counter % 10000 == 0:
                print(f"{counter} of {n} :::: {counter/n*100:5.2f}% ::: {len(result_As)} finded combinations")
            As = comb[::2]
            e = comb[1::2]
            
            Qav = list(map(lambda x: np.sum(self.Qav_As[:] * As) + np.sum(self.Qav_e[:] * e) + x.Qav, self.cyclograms))
            T = list(map(lambda x: (x/(sigma*np.sum(self.areas[:]*e)))**0.25, Qav))
            if sum(list(map(lambda x: x < Tmax and x > Tmin, T))) == len(T):
                result_As.append(As)
                result_e.append(e)
            counter += 1
        
        df = pd.concat([pd.DataFrame(result_As), pd.DataFrame(result_e)], axis=1)
        new_columns_name_As = [f"As_{i}" for i in range(6)]
        new_columns_name_e = [f"e_{i}" for i in range(6)]
        new_columns_name = new_columns_name_As + new_columns_name_e
        df.columns = new_columns_name
        
        return df, new_columns_name_As, new_columns_name_e        