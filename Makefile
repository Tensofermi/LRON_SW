# 设置最终目标文件名
TARGET := a.out

# 定义目录
SRC_DIR := src

# SRC_DIR := src_2D_Ising_v1

# SRC_DIR := src_2D_Ising_v2

BIN_DIR := bin
BUILD_DIR := build

# 查找所有源代码文件（.cpp）
SRC := $(shell find $(SRC_DIR) -type f -name *.cpp)

# 定义目标文件路径，源文件 (.cpp) 转换为目标文件 (.o)
OBJ := $(SRC:$(SRC_DIR)/%.cpp=$(BUILD_DIR)/%.o)

# 定义依赖文件路径 (.d)
DEP := $(OBJ:.o=.d)

# 查找所有头文件目录
INC_DIR := $(shell find $(SRC_DIR) -type d)

# 生成包含所有头文件目录的 -I 标志
INC_FLAG := $(addprefix -I,$(INC_DIR))

# 设置编译器
CXX := g++

# 定义所需的库
ARMADILLO_LIB := -larmadillo
LAPACK_LIB := -llapack
BLAS_LIB := -lblas
GMP_LIB := -lgmp
GMPXX_LIB := -lgmpxx
GSL_LIB := -lgsl
BOOST_LIB := -lboost_mpi -lboost_serialization -lboost_filesystem -lboost_system
FFT_LIB := 

# 定义最终链接的所有库
ALL_LIBS := $(FFT_LIB)

# 编译器警告和调试选项
CXX_WARNINGS := -Wall -pedantic-errors
CXX_DEBUG :=  # 可以在这里启用调试选项
CXX_OPTIM := -O3  # 优化等级
CXX_DEP := -MMD  # 生成依赖文件


CXXFLAGS := -std=c++11 $(CXX_DEBUG) $(CXX_WARNINGS) $(CXX_OPTIM) $(CXX_DEP) $(INC_FLAG)


$(BIN_DIR)/$(TARGET): $(OBJ)
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $(OBJ) -o $@ $(ALL_LIBS)


$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) -c -o $@ $<  


.PHONY: clean
clean:
	$(RM) -rf $(BUILD_DIR)/* $(BIN_DIR)/$(TARGET) $(BUILD_DIR) $(BIN_DIR)  # 删除构建目录和可执行文件


-include $(DEP)
