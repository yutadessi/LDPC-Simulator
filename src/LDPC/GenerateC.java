package LDPC;

import java.util.Random;
public class GenerateC {
    public static int[] geneC(int [][] encodedG){
        //元のメッセージkと送信語cの初期化
        int[] k = new int[encodedG.length];
        int[] c = new int[encodedG[0].length];
        Random random = new Random();

        //元のメッセージ作成
        for(int i = 0;i < k.length;i++){
            k[i] = random.nextInt(2);
        }
//        System.out.print("メッセージk ");
//        Print.Array(k);

        //送信語cの作成
        for(int i = 0;i < encodedG[0].length;i++){
            for(int j = 0;j < k.length;j++){
                c[i] += k[j] * encodedG[j][i];
            }
        }
        for(int i = 0;i < c.length;i++){
            c[i] %= 2;
        }
//        System.out.print("送信語c ");
//        Print.Array(c);

        return c;
    }
}
