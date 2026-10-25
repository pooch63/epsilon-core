import type { Metadata } from "next";
import "./globals.css";

export const metadata: Metadata = {
  title: "Epsilon-core viewer",
  description: "Look up heuristic node numbers in an Amazon reviews JSONL and show the user–product edges.",
};

export default function RootLayout({ children }: { children: React.ReactNode }) {
  return (
    <html lang="en">
      <body>{children}</body>
    </html>
  );
}
