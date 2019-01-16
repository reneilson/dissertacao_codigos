def classifier(soma_fft, max_, rms_dados, 
               std_, dist_max_mean_dados,
               dist_min_mean_dados, dist_cont,
               min_env_dados, max_env_dados,
               dist_env_fft_mx, dist_env_fft_mx_mean,
               max_aux):
    if ( min_env_dados <= 0.9993405342102051 ) :
        if ( rms_dados <= 0.688856840133667 ) :
            if ( min_env_dados <= 0.10013996809720993 ) :
                return 1
            else:
                if ( dist_env_fft_mx <= 2.1340372562408447 ) :
                    if ( dist_min_mean_dados <= 1.073239803314209 ) :
                        if ( std_ <= 0.2932622730731964 ) :
                            return 2
                        else:
                            return 1
                        
                    else:
                        if ( dist_env_fft_mx_mean <= 0.7787008881568909 ) :
                            return 2
                        else:
                            return 1
                        
                    
                else:
                    if ( max_**2 <= 18764.609375 ) :
                        return 1
                    else:
                        if ( soma_fft <= 59.26722717285156 ) :
                            return 2
                        else:
                            return 1
                        
                    
                
            
        else:
            if ( max_**2 <= 0.00017393683083355427 ) :
                if ( soma_fft <= 4.102768898010254 ) :
                    if ( dist_max_mean_dados <= 0.9998793005943298 ) :
                        if ( max_env_dados <= 0.9593193531036377 ) :
                            return 6
                        else:
                            return 0
                        
                    else:
                        return 5
                    
                else:
                    if ( max_env_dados <= 1.0228557586669922 ) :
                        if ( max_env_dados <= 0.971593976020813 ) :
                            return 6
                        else:
                            return 0
                        
                    else:
                        return 6
                    
                
            else:
                if ( max_env_dados <= 1.0003225803375244 ) :
                    if ( soma_fft <= 7.246959686279297 ) :
                        return 2
                    else:
                        return 8
                    
                else:
                    if ( dist_env_fft_mx_mean <= 1.623831033706665 ) :
                        return 9
                    else:
                        return 6
                    
                
            
        
    else:
        if ( max_aux <= 2.5 ) :
            if ( dist_min_mean_dados <= 1.019820213317871 ) :
                if ( max_**2 <= 2.6601304625728517e-07 ) :
                    if ( dist_env_fft_mx_mean <= 0.0754564106464386 ) :
                        if ( dist_max_mean_dados <= 1.0000007152557373 ) :
                            return 8
                        else:
                            return 5
                        
                    else:
                        if ( dist_min_mean_dados <= 1.0012495517730713 ) :
                            return 0
                        else:
                            return 8
                        
                    
                else:
                    if ( dist_min_mean_dados <= 1.000420331954956 ) :
                        if ( dist_min_mean_dados <= 1.0000941753387451 ) :
                            return 4
                        else:
                            return 7
                        
                    else:
                        if ( soma_fft <= 6.391940116882324 ) :
                            return 3
                        else:
                            return 8
                        
                    
                
            else:
                if ( rms_dados <= 0.720308780670166 ) :
                    if ( rms_dados <= 0.7111003398895264 ) :
                        return 4
                    else:
                        if ( max_**2 <= 5818.8583984375 ) :
                            return 6
                        else:
                            return 9
                        
                    
                else:
                    return 3
                
            
        else:
            if ( min_env_dados <= 1.0003288984298706 ) :
                if ( dist_env_fft_mx_mean <= 0.05461488664150238 ) :
                    if ( dist_min_mean_dados <= 1.0001726150512695 ) :
                        return 4
                    else:
                        if ( dist_min_mean_dados <= 1.0079435110092163 ) :
                            return 7
                        else:
                            return 4
                        
                    
                else:
                    return 4
                
            else:
                return 6