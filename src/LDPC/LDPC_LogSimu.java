package LDPC;

import java.util.Arrays;
import java.util.List;
import java.util.ArrayList;
import java.io.PrintWriter;
import java.io.IOException;
import java.nio.charset.Charset;

//--------------------デスクトップPCでの処理速度--------------------
//処理速度(フレーム数:10000,1024-8-4サイズ,誤り率数:10,Lmax:20 ):18分
//処理速度(フレーム数:10000,1024-8-4サイズ,誤り率数:10,Lmax:100):31分

public class LDPC_LogSimu {
    public static void main(String[] args) {

        //ファイル名、毎回変える！！--------
        String fileNAMEME = "New8-(3-5)";
        //------------------------------
        String fileNames = fileNAMEME + "-result.csv"; //結果保存ファイル名
        String filePath = fileNAMEME + "-HMatrix.txt"; //検査行列保存ファイル名

        //符号パラメータ
        int n = 1024; //符号長
        int wr = 8; //行重み(n % wr = 0)
        int[] wcs = {3,4,5}; //列重み
        int maxL = 50; //最大反復回数
        int numFrames = 10000; //フレーム数

        //通信路誤り率eの集合
        double[] eValues = {0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1};
//            double[] eValues = {0.05};

        //出力用保存配列
        //-符号パラメータ
        int[][] CodeParameter = new int[wcs.length][5];
        //-各誤り率のデータ
        double[][][] channelBitErrorRate = new double[wcs.length][eValues.length][numFrames]; //各フレームの実際の誤り率
        double[][][] frameErrorRate = new double[wcs.length][eValues.length][numFrames]; //各フレームのFER
        double[][][] averageTrueIterations = new double[wcs.length][eValues.length][numFrames]; //訂正成功時の平均繰り返し回数
        double[][][] averageFalseIterations = new double[wcs.length][eValues.length][numFrames]; //訂正失敗時の平均繰り返し回数
        double[][] residualsErrorBits = new double[wcs.length][eValues.length]; //情報ビットの残留誤りビット数
        double[][] errorCollectionBits = new double[wcs.length][eValues.length]; //情報ビットの誤訂正ビット数
        double[][][] iterationDistribution = new double[wcs.length][eValues.length][maxL]; //反復回数の度数分布

        for(int wc : wcs){

            //検査行列Hと生成行列Gの作成
            int [][] H = GenerateMatrix.gallagerCheckMatrix(n,wr,wc);
            int [][] G = GenerateMatrix.generatorMatrix(H);

            //検査行列を保存
            CheckMatrixIO.saveCheckMatrix(H,filePath);

            //HとGを組織符号化
            List<Integer> columnIndicatesToSwap = new ArrayList<>();  //組織符号化用インデックス
            int[][] encodedG = TissueEncoder.EncodeG(G,columnIndicatesToSwap);
            int[][] encodedH = TissueEncoder.EncodeH(H,columnIndicatesToSwap);

            //復号
            for(double e : eValues){

                for(int frame = 0;frame < numFrames;frame++){
                    //メッセージと送信語、受信語の作成
                    int[] c = GenerateC.geneC(encodedG);
                    int[] r = Channel.GenerateR(c,e);

                    //実際の通信路での誤り率の取得
                    double cBER = Channel.CheckError(c,r);

                    //対数領域sum-product復号&Min-Sum復号法
                    LogDecoder.DecodeResult result = LogDecoder.decode(encodedH,r,e,maxL);
//                    MinSumDecoder.DecodeResult result = MinSumDecoder.decode(encodedH,r,e,maxL);

                    int[] decodedC = result.decodedCode();
                    int iterations = result.iterationNum();

                }
            }
        }

        try (PrintWriter pw = new PrintWriter(fileNames, Charset.forName("Windows-31j"))){
//            pw.printf("%s,%s,%s,%s,%s\n","符号長","行重み","列重み","最大反復回数","フレーム数");
//            pw.printf("%s,%s,%s,%s,%s\n\n",n,wr,wc,maxL,numFrames);
//            pw.printf("%s,%s,%s,%s,%s,%s,%s\n","通信路誤り率","実際の通信路誤り率の平均","FER","IBER","成功時の平均繰り返し回数","失敗時の平均繰り返し回数","平均誤訂正ビット数");
//            pw.printf("%.2f,%s,%.4f,%s,%s,%s,%s\n", e, aveCBER, fer, iber,aveTrueIterations,aveFalseIterations,aveMissCorrection);
//            pw.printf("\n%s,%s,%s,%s\n","C-BER","Frame","ErrorIBits","Iterations");
//            pw.printf("%s,%s,%s,%s\n",groupOfCBER[i][j],groupOfFrame[i][j],groupOfErrorInfoBits[i][j],groupOfIterations[i][j]);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
