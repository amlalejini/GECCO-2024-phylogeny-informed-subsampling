# SignalGP instruction set

Below, we document the instruction set used in our GP system for our 2023 GPTP experiments.

Abbreviations:

- EOP: End of program
- Reg: local register
  - Reg[0] indicates the value at the register specified by an instruction's first _argument_, Reg[1] indicates the value at the register specified by an instruction's second argument, and Reg[2] indicates the value at the register specified by the instruction's third argument.
  - Reg[0], Reg[1], _etc_: Register 0, Register 1, _etc._
- Input: input buffer
  - Follows same scheme as Reg
- Output: output buffer
  - Follows same scheme as Reg
- Global: global memory buffer
  - Follows same scheme as Reg
- Arg: Instruction argument
  - Arg[i] indicates the i'th instruction argument (an integer encoded in the genome)
  - E.g., Arg[0] is an instruction's first argument

Instructions that would produce undefined behavior (e.g., division by zero) are treated as no operations.

## Default Instructions

I.e., instructions used across all diagnostic tasks.


| Instruction | Arguments Used | Description |
| :--- | :---: | :--- |
| `Nop` | 0 | No operation |
| `Not` | 1 | Reg[0] = !Reg[0] |
| `Inc` | 1 | Reg[0] = Reg[0] + 1 |
| `Dec` | 1 | Reg[0] = Reg[0] - 1 |
| `Add` | 3 | Reg[2] = Reg[0] + Reg[1] |
| `Sub` | 3 | Reg[2] = Reg[0] - Reg[1] |
| `Mult`  | 3 | Reg[2] = Reg[0] * Reg[1] |
| `Div` | 3 | Reg[2] = Reg[0] / Reg[1] |
| `Mod` | 3 | Reg[2] = Reg[0] % Reg[1] |
| `Nand` | 2 | Reg[2] = !(Reg[0] & Reg[1]) |
| `TestEqu`  | 3 | Reg[2] = Reg[0] == Reg[1] |
| `TestNEqu` | 3 | Reg[2] = Reg[0] != Reg[1] |
| `TestLess` | 3 | Reg[2] = Reg[0] < Reg[1] |
| `TestLessEqu` | 3 | Reg[2] = Reg[0] <= Reg[1] |
| `TestGreater`  | 3 | Reg[2] = Reg[0] > Reg[1] |
| `TestGreaterEqu`  | 3 | Reg[2] = Reg[0] >= Reg[1] |
| `SetMem` | 2 | Reg[0] = Arg[1] |
| `Terminal`  | 1 | Reg[0] = double value encoded by instruction tag |
| `CopyMem` | 2 | Reg[0] = Reg[1] |
| `SwapMem` | 2   | Swap(Reg[0], Reg[1]) |
| `InputToWorking`  | 2    | Reg[1] = Input[0] |
| `WorkingToOutput`  | 2    | Output[1] = Reg[0] |
| `If`      | 1   | If Reg[0] != 0, proceed. Otherwise skip to the next `Close` or EOP. |
| `While` | 1     | While Reg[0] != 0, loop. Otherwise skip to next `Close` or EOP. |
| `Close` | 0     | Indicate the end of a control block of code (e.g., loop, if). |
| `Break` | 0     | Break out of current control flow (e.g., loop). |
| `Call`  | 0    | Call a function, using this instruction's tag to determine which function is called. |
| `Routine`  | 0  | Same as call, but local memory is shared. Sort of like a jump that will jump back when the routine ends. |
| `Return`  | 0    | Return from the current function call. |
| `WorkingToGlobal` | 2 | Global[1] = Reg[0] |
| `GlobalToWorking` | 2 | Reg[1] = Global[0] |
| `FullGlobalToWorking` | 0 | Copy entire global memory buffer into working memory buffer |
| `FullWorkingToGlobal` | 0 | Copy entire working memory buffer into global memory buffer |

Note that `Nand` performs a bitwise operation.

## Problem-specific instructions

Each problem has problem-specific instructions for producing output.

### Fizz Buzz

- SubmitFizz
- SubmitBuzz
- SubmitFizzBuzz
- SubmitEcho

### Median

- SubmitOutput

### Grade

- SubmitA
- SubmitB
- SubmitC
- SubmitD
- SubmitF

### Small or large

- SubmitSmall
- SubmitLarge
- SubmitNeither