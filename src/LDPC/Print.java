package LDPC;

import java.io.PrintWriter;

public class Print {
    public static void Matrix(int[][] Mat, PrintWriter pw){
        pw.println("Matrix [" + Mat.length + ", " + Mat[0].length + "]");
        for(int i = 0;i < Mat.length;i++){
            for(int j = 0;j < Mat[0].length;j++) {
                pw.print(Mat[i][j] + " ");
            }
            pw.println(" ");
        }
        pw.println(" ");
    }

    public static void DetailMatrix(int[][] Mat, PrintWriter pw){
        pw.print("encodedH = {");
        for(int i = 0;i < Mat.length;i++){
            pw.print("{");
            for(int j = 0;j < Mat[0].length;j++) {
                if(j != Mat[0].length-1){
                    pw.print(Mat[i][j] + ",");
                }else{
                    pw.print(Mat[i][j]);
                }
            }
            if(i != Mat.length-1){
                pw.print("},");
            }else{
                pw.print("}");
            }
        }
        pw.println("}");
    }

    public static void Array(int[] Array,PrintWriter pw){
        pw.println("Array [" + Array.length + "]");
        for(int i = 0;i < Array.length;i++){
            pw.print(Array[i] + " ");
        }
        pw.println(" \n");
    }
}
