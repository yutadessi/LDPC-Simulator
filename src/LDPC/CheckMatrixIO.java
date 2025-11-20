package LDPC;

import java.util.Scanner;
import java.io.*;

public class CheckMatrixIO {
    //検査行列を保存するメソッド
    public static void saveCheckMatrix(int[][] H,String filePath){
        try (PrintWriter pw = new PrintWriter(new FileWriter(filePath))) {
            int m = H.length;
            int n = H[0].length;

            pw.println(m + " " + n);
            for(int i = 0;i < m;i++){
                for(int j = 0;j < n;j++){
                    if(H[i][j] == 1){
                        pw.print(i + " " + j + " ");
                    }
                }
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    //検査行列を読み取り、復元するメソッド
    public static int[][] loadCheckMatrix(String filePath){
        try(Scanner scanner = new Scanner(new File(filePath))){

            if(!scanner.hasNextInt()) return null;

            int m = scanner.nextInt();
            int n = scanner.nextInt();

            int[][] H = new int[m][n];

            while(scanner.hasNextInt()){
                int r = scanner.nextInt();
                if (scanner.hasNextInt()) {
                    int c = scanner.nextInt();
                    if(r < m && c < n){
                        H[r][c] = 1;
                    }
                }
            }

            return H;
        } catch (FileNotFoundException e) {
            e.printStackTrace();
            return null;
        }
    }
}
