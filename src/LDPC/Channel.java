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
    public static void CheckError (int[] c,int[] r,double e){
        double count = 0;
        for(int i = 0;i < c.length;i++){
            count += c[i] - r[i] != 0 ? 1 : 0;
        }
        System.out.println("通信路誤り率e = " + e + ",実際の通信路誤り率:" + (count / c.length));
    }
}
