package LDPC;

import java.util.List;
import java.util.ArrayList;
public class LDPC_SimuGT {
    public static void main(String[] args) {
        //n:符号長、wr:行重み、wc:列重みの初期化
        int n = 18;
        int wr = 6;
        int wc = 3;

        //通信路誤り率eの設定
        double e = 0.1;

        //最大反復回数maxLを設定
        int maxL = 20;

        //組織符号化する際の移動列インデックス
        List<Integer> columnIndicatesToSwap = new ArrayList<>();

        //行重みが符号長に対して適切か確認
        if(Divisibility.checker(n,wr) == 1) return;

        //検査行列H作成と表示
        int [][] H = GenerateMatrix.gallagerCheckMatrix(n,wr,wc); //検査行列H
        System.out.print("H ");
        Print.Matrix(H);

        //生成行列G作成と表示
        int [][] G = GenerateMatrix.generatorMatrix(H,n,wr,wc);
        System.out.print("G ");
        Print.Matrix(G);

        //適切な生成行列か確認
        ConfirmG.ConG(G,H);

        //組織符号化生成とその確認
        int[][] encodedG = TissueEncoder.EncodeG(G,columnIndicatesToSwap);
        int[][] encodedH = TissueEncoder.EncodeH(H,columnIndicatesToSwap);
        System.out.println(columnIndicatesToSwap+"\n");
        System.out.print("Encoded G ");
        Print.Matrix(encodedG);
        System.out.print("Encoded H ");
        Print.Matrix(encodedH);
        ConfirmG.ConG(encodedG,encodedH);

        //メッセージと送信語の作成と表示
        int[] c = GenerateC.geneC(encodedG);

        //受信語作成
        int[] r = Channel.GenerateR(c,e);
        Channel.CheckError(c,r,e);

        //確率領域sum-product復号
        int[] probabilityDecodedWords = ProbDecoder.decode(encodedH,r,e,maxL);


    }
}
