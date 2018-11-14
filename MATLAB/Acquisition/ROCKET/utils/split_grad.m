function   [Gxs_rf, Gys_rf, Gzs_rf, Gxs_adc,Gys_adc,Gzs_adc] = split_grad(Gxs,Gys,Gzs,Grf)
            %% RF bit
             Gxs_rf = Gxs;
             Gxs_rf.waveform = Gxs_rf.waveform(1:Grf);
             Gxs_rf.t = Gxs_rf.t(1:Grf);
             
             Gys_rf = Gys;
             Gys_rf.waveform = Gys_rf.waveform(1:Grf);
             Gys_rf.t = Gys_rf.t(1:Grf);
             
             Gzs_rf = Gzs;
             Gzs_rf.waveform = Gzs_rf.waveform(1:Grf);
             Gzs_rf.t = Gzs_rf.t(1:Grf);
             
             %% ADC bit
             Gxs_adc = Gxs;
             Gxs_adc.waveform = Gxs_adc.waveform(Grf+1:end);
%              Gxs_adc.t = Gxs_adc.t(Grf+1:end);
              Gxs_adc.t = Gxs_adc.t(Grf+1:end);
             
             Gys_adc = Gys;
             Gys_adc.waveform = Gys_adc.waveform(Grf+1:end);
             Gys_adc.t = Gys_adc.t(Grf+1:end);
             
             Gzs_adc = Gzs;
             Gzs_adc.waveform = Gzs_adc.waveform(Grf+1:end);
             Gzs_adc.t = Gzs_adc.t(Grf+1:end);