      if (method==1)
          time_FE;
      elseif (method==2)
          time_AB_CN;
      elseif (method==31)
          time_BE;
      elseif (method==32)
          time_BE_u;
      elseif (method==4)
          time_CN;
      elseif (method==5)
          time_oneleg;
      elseif (method==61)
%           if (use_Schur==0)
% %             time_FC1;
%             time_IRK2_linear;
%           elseif (use_Schur==1)
%             time_FC1_Schur;
%           end
%       elseif (method==62)
%           if (use_Schur==0)
% %               time_FC2;
%               time_IRK2_u;
%           elseif (use_Schur==1)
%               time_FC2_Schur;
%           end
        if (order4==0)
          time_Gauss2_u;
        elseif (order4==1)
          time_Gauss2_u_4thorder;
        end
      elseif (method==62)
          time_Gauss2_u_linear;
      elseif (method==63)
          time_IRK2_u_PC;          
      elseif (method==71)
          time_IM1;
      elseif (method==72)
          time_IM2;
      elseif (method==81 || method==82)
          time_RK4;
      elseif (method==91)
          if (order4==0)
            time_Gauss4_u;
          elseif (order4==1)
            time_Gauss4_u_4thorder;
          end
      elseif (method==92)
          time_Gauss4_u_linear;
%       elseif (method==93)
%           time_IRK4_u_PC;
%           time_IRK4_p;
%           time_IRK4_PC;
      elseif (method==101 || method==102)
          time_RK2;
      elseif (method==111)
          time_RK3_CN;
      elseif (method==112)
          time_RK3;
      elseif (method==12)
          time_RK_gen3;
%           time_RK_gen3_movingAD;
%           time_RK_gen4;
%           time_RK_diffusion;
      elseif (method==13)
          time_SDIRK;
      elseif (method==141)
          time_ARK_RadauIIAB;    
      elseif (method==142)
          time_ARK_RadauIIAB_linear;
      elseif (method==15)
          time_Gauss6_u;
      elseif (method==16)
          time_LobIIIC_u;
      elseif (method==171)
          time_ARK_DIRK;          
      elseif (method==172)
          time_ARK_DIRK_linear;          
      elseif (method==181)
          time_ARK_LobIIICE;
      elseif (method==182)
          time_ARK_LobIIICE_linear; 
      elseif (method==191)
          time_CN_u; % 2-stage
%           time_LobIIIA_u; % 3-stage
%           time_RadIIA_u_s3;
      elseif (method==192)
          time_CN_u_linear;
      end