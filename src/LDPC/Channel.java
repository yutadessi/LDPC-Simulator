package LDPC;

import java.util.Random;
public class Channel {
    public static int[] GenerateR (int[] c,double e){
        int[] r = new int[c.length];
        Random random = new Random();
        for(int i = 0;i < r.length;i++){
            int rand = random.nextInt(1000);
            if(rand < e * 1000){
                r[i] = c[i] == 1 ? 0 : 1;
            }else {
                r[i] = c[i];
            }
        }
        return r;
    }
    public static double CheckError(int[] c,int[] r){
        int countCBE = 0;
        for(int i = 0;i < c.length;i++){
            if(c[i] != r[i]){
                countCBE ++;
            }
        }
        double cBER = (double) countCBE/c.length;
        return cBER;
    }

}
