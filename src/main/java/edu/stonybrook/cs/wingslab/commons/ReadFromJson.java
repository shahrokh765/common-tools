package edu.stonybrook.cs.wingslab.commons;

import java.io.*;
import java.util.Map;
import java.lang.reflect.Type;
import java.util.concurrent.ConcurrentHashMap;

import com.google.gson.Gson;
import com.google.gson.reflect.TypeToken;

public class ReadFromJson {
    public static void main(String args[]) throws FileNotFoundException {
        String tt = "{'10':15.2, '12':14.4}";
        Type mapOfStringObjectType = new TypeToken<ConcurrentHashMap<String,ConcurrentHashMap<String, Double>>>() {}.getType();
        Gson gson = new Gson();
//        ConcurrentHashMap<String, Double> obj1 = gson.fromJson(tt, mapOfStringObjectType);
//        System.out.println(System.getProperty("user.dir"));
//        new FileReader("C:\\fileName.json")
        ConcurrentHashMap<String, ConcurrentHashMap<String, Double>> obj =
                gson.fromJson(new FileReader("resources/splat/pl_map/pl_map.json"), mapOfStringObjectType);
//        System.out.println(obj);
        try(FileOutputStream fos = new FileOutputStream("resources/splat/pl_map/pl_map.sr");
            ObjectOutputStream oos = new ObjectOutputStream(fos)) {
            oos.writeObject(obj);
        } catch (IOException e) {
            e.printStackTrace();
            System.out.println("Serialization was not successful! Saving map failed");
        }
    }
}