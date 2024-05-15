dir_saving = getDirectory("Choose a Directory to save");

dir_processing = getDirectory("Choose a Directory to proess"); 

list = getFileList(dir_processing);

for (i = 0; i < list.length; i++) {

  open("D:/360MoveData/Users/Haofei Ge/Desktop/该读书了/细胞周期伪时间课题（毕业论文版）/活细胞与伪时间的验证数据/伪时间处理/preprocess/Input/" + list[i]);
  
  run("Subtract Background...", "rolling=50"); 
  
  saveAs("bmp", dir_saving + list[i]);

  close();
}
