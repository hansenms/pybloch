import numpy as np
import string
import cmath

def read_poet_simulation(filename):
    header_lines = 0
    with open(filename,'r') as f:
        for line in f:
            header_lines = header_lines + 1
            l = string.strip(line)
            if len(l) == 0:
                continue

            if l[0] == ';':
                if string.find(line, 'ADC Signal Data') > 0:
                    #This is where we read the field names
                    field_names = line[1:].split("\t")
                    break
                else:
                    continue

        field_names = map(string.strip, field_names)
        f.close()

        sequence_table_raw = np.loadtxt(filename, skiprows=header_lines)

        sequence_table = np.zeros((sequence_table_raw.shape[0],5),dtype=np.complex64)

        #Now convert the table
        rf_signal_column = field_names.index('RF-Signal (ch. 0, 1H, 63.6 MHz)')
        rf_phase_column = field_names.index('RF-Signal Phase (ch. 0, 1H, 63.6 MHz)')
        nco_phase_column = field_names.index('Numeric Crystal Oscillator 1 Phase')
        z_gradient_column = field_names.index('Z Gradient (GPA 0)')
        y_gradient_column = field_names.index('Y Gradient (GPA 0)')
        x_gradient_column = field_names.index('X Gradient (GPA 0)')
        adc_column = field_names.index('ADC Signal Data')

        #Read values and convert to SI units
        for t in range(sequence_table_raw.shape[0]):
            if (sequence_table_raw[t,rf_signal_column] > 0):
                #TODO: Have to figure out correct RF scaling
                sequence_table[t,0] = 0.037570*1.0e-6*cmath.rect(sequence_table_raw[t,rf_signal_column], np.pi*((sequence_table_raw[t,rf_phase_column] + sequence_table_raw[t,nco_phase_column])/180.0))
            sequence_table[t,1] = 1.0e-3*sequence_table_raw[t,z_gradient_column]
            sequence_table[t,2] = 1.0e-3*sequence_table_raw[t,y_gradient_column]
            sequence_table[t,3] = 1.0e-3*sequence_table_raw[t,x_gradient_column]
            sequence_table[t,4] = cmath.rect(sequence_table_raw[t,adc_column], np.pi*(sequence_table_raw[t,nco_phase_column]/180.0))

        return sequence_table
