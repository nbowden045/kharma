if [[ "$OSTYPE" == "darwin"* ]]; then
    BREW_LLVM="$(brew --prefix llvm)"
    export LDFLAGS="$LDFLAGS -L$BREW_LLVM/lib/c++"
    C_NATIVE="$BREW_LLVM/bin/clang"
    CXX_NATIVE="$BREW_LLVM/bin/clang++"
fi
