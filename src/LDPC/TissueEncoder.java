package LDPC;
import java.util.List;
import java.util.ArrayList;

public class TissueEncoder {
    //生成行列の組織符号化
    public static int[][] EncodeG(int[][] G,List<Integer> column){
        //生成行列Gのコピー作成
        int[][] enG = new int[G.length][G[0].length];
        for(int i = 0;i < G.length;i++)enG[i] = G[i].clone();

        //列ごとにガウスの消去法を適用
        for(int i = 0;i < enG.length;i++){
            int[] tempC = new int[enG.length];
            int[] tempR;
            List<Integer> Row = new ArrayList<>();
            List<Integer> Col = new ArrayList<>();
            //列交換処理
            if(enG[i][i] != 1){
                for(int jc = 0;jc < enG[0].length;jc++){
                    if(enG[i][jc] == 1){
                        Col.add(jc);
                        column.add(i);
                        column.add(jc);
                        break;
                    }
                }
                for(int jc2 = 0;jc2 < enG.length;jc2++){
                    tempC[jc2] = enG[jc2][Col.get(0)];
                    enG[jc2][Col.get(0)] = enG[jc2][i];
                    enG[jc2][i] = tempC[jc2];
                }
                Col.clear();
            }
            //加算処理
            for(int jr = 0;jr < enG.length;jr++){
                if(enG[jr][i] == 1 && jr != i)Row.add(jr);
            }
            if(Row.size() != 0){
                for(int jr2 = 0;jr2 < Row.size();jr2++){
                    for(int addRow = 0;addRow < enG[0].length;addRow++){
                        if(enG[i][addRow] == 1 && enG[Row.get(jr2)][addRow] == 1){
                            enG[Row.get(jr2)][addRow] = 0;
                            continue;
                        }
                        if(enG[i][addRow] == 1 && enG[Row.get(jr2)][addRow] == 0)enG[Row.get(jr2)][addRow] = 1;
                    }
                }
            }
            Row.clear();
        }
        return enG;
    }

    //検査行列の組織符号化
    public static int[][] EncodeH(int[][] H,List<Integer> column){
        //検査行列のコピー作成
        int[][] enH = new int[H.length][H[0].length];
        for(int i = 0;i < H.length;i++)enH[i] = H[i].clone();

        //検査行列の列交換
        int[] temp = new int[enH.length];
        for(int i = 0;i < (column.size() / 2);i++){
            for(int j = 0;j < enH.length;j++){
                temp[j] = enH[j][column.get(2 * i)];
                enH[j][column.get(2 * i)] = enH[j][column.get((2 * i) + 1)];
                enH[j][column.get((2 * i) + 1)] = temp[j];
            }
        }
        return enH;
    }
}
