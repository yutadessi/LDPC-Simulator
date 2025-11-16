package LDPC;

import java.io.PrintWriter;

public class ConfirmG {
    public static void ConG(int[][] G,int[][] H,PrintWriter pw){
        int total[][] = new int[G.length][H.length];
        int parityCheck = 0;
        for(int i = 0;i < G.length;i++){
            for(int j = 0;j < H.length;j++){
                for(int k = 0;k < H[0].length;k++){
                    total[i][j] += G[i][k] * H[j][k];
                }
            }
        }
        System.out.print("G * Ht ");
        Print.Matrix(total,pw);
        for(int[] i  : total){
            for(int j : i){
                parityCheck += j % 2;
            }
        }
        if(parityCheck == 0) System.out.println("お使いの生成行列は正常です。\n");
        if(parityCheck != 0) System.out.println("お使いの生成行列は異常です。\n");
    }
}
