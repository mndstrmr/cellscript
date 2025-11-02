module hello_world(
  input logic clk_i,
  output logic exit_o,
  output logic display_o,
  output logic [3:0] char_o
);

logic [3:0] time_q;

always @(posedge clk_i) time_q += 1;

assign exit_o = time_q == 4;
assign display_o = 1;

always_comb case (time_q)
  0: char_o = 7;
  1: char_o = 4;
  2: char_o = 11;
  3: char_o = 11;
  4: char_o = 14;
  default: char_o = 0;
endcase

endmodule
