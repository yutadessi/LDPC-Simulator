package LDPC;

import java.util.Arrays;
import java.util.List;
import java.util.ArrayList;
public class LDPC_SimuGT {
    public static void main(String[] args) {

        //4fps
        //符号パラメーター
        int n = 1024; //符号長
        int wr = 8; //行重み
        int wc = 4; //列重み
        int maxL = 20; //最大反復回数

        //シミュレーション設定
        int numFrames = 500;

        //通信路誤り率eの設定
        double[] eValues = {0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1};

        //行重みが符号長に対して適切か確認
//        if(Divisibility.checker(n,wr) == 1) return;

        //検査行列Hと生成行列Gの作成と表示
        int [][] H = GenerateMatrix.gallagerCheckMatrix(n,wr,wc);
//        System.out.print("H ");
//        Print.Matrix(H);
        int [][] G = GenerateMatrix.generatorMatrix(H,n,wr,wc);
//        System.out.print("G ");
//        Print.Matrix(G);
//        ConfirmG.ConG(G,H); //適切な生成行列か確認


        //HとGを組織符号化
        List<Integer> columnIndicatesToSwap = new ArrayList<>();
        int[][] encodedG = TissueEncoder.EncodeG(G,columnIndicatesToSwap);
        int[][] encodedH = TissueEncoder.EncodeH(H,columnIndicatesToSwap);
//        System.out.println(columnIndicatesToSwap+"\n");
//        System.out.print("Encoded G ");
//        Print.Matrix(encodedG);
//        System.out.print("Encoded H ");
//        Print.Matrix(encodedH);
//        ConfirmG.ConG(encodedG,encodedH); //適切な生成行列か確認

        for(double e : eValues){
            long frameErrorCount = 0;

            for(int frame = 0;frame < numFrames;frame++){

                //メッセージと送信語の作成と表示
                int[] c = GenerateC.geneC(encodedG);

                //受信語作成
                int[] r = Channel.GenerateR(c,e);

                //確率領域sum-product復号
                int[] estimatedC = ProbDecoder.decode(encodedH,r,e,maxL);

                if(!Arrays.equals(c,estimatedC)){
                    frameErrorCount++;
                }
            }
            double fer = (double)frameErrorCount/numFrames;
            System.out.printf("%.2f \t\t| %.6f\n",e,fer);
        }


    }
}
