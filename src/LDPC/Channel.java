package LDPC;

import java.util.Random;
public class Channel {
    public static int[] GenerateR (int[] c,double e,int gLength){
        int[] r = new int[c.length];
        Random random = new Random();
        for(int i = 0;i < gLength;i++){
            if(random.nextDouble() < e){
                r[i] = c[i] == 1 ? 0 : 1;
            }else {
                r[i] = c[i];
            }
        }
        for(int i = gLength;i < r.length;i++)r[i] = c[i];
        return r;
    }
    public static double CheckError(int[] c,int[] r,int gLength){
        int countCBE = 0;
        for(int i = 0;i < gLength;i++){
            if(c[i] != r[i]){
                countCBE ++;
            }
        }
        double cBER = (double) countCBE/gLength;
        return cBER;
    }
    public static int[] GenerateRBefore (int[] c,double e,int gLength){
        int[] r = new int[c.length];
        Random random = new Random();
        for(int i = 0;i < r.length;i++){
            if(random.nextDouble() < e){
                r[i] = c[i] == 1 ? 0 : 1;
            }else {
                r[i] = c[i];
            }
        }
        return r;
    }

}
