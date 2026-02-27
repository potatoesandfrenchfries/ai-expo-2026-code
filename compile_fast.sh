#!/bin/bash
# Fast compilation script for optimized molecular modules

echo "=== Fast Compilation for Molecular Chemistry ==="
echo ""

# Regenerate Makefile from _CoqProject
echo "Regenerating Makefile..."
coq_makefile -f _CoqProject -o Makefile

echo ""
echo "=== Compiling optimized modules only ==="
echo ""

# Compile only the fast modules
echo "1. Compiling Geometry..."
rocq compile -q src/rocq/Geometry.v

echo "2. Compiling Bonds..."
rocq compile -q src/rocq/Bonds.v

echo "3. Compiling CommonElements (10 elements instead of 118)..."
rocq compile -q src/rocq/CommonElements.v

echo "4. Compiling FastMolecule..."
rocq compile -q src/rocq/FastMolecule.v

echo "5. Compiling FastDemo (ethanol + aspirin)..."
rocq compile -q src/rocq/FastDemo.v

echo ""
echo "=== Compilation complete! ==="
echo ""
echo "To test the molecules, run:"
echo "  rocq repl -R src/rocq Chemistry"
echo ""
echo "Then in the REPL:"
echo "  Require Import Chemistry.FastDemo."
echo "  Eval native_compute in (is_connected aspirin)."
echo "  Eval native_compute in (molecular_weight aspirin)."
echo ""
