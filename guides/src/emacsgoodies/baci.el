;; baci.el
;; $Id: baci.el,v 1.1 2006/11/20 01:16:04 pbb Exp $
;;
;; About: A simple Emacs mode for Baci input DAT files
;;
;; Typist: B Bornemann
;; Notice: File comes without any warranty what so ever.
;;
;; Based on: Jorgen Bergstrom's lsdyna.el <jorgen@polymerFEM.com>
;; URL:      http://polymerFEM.com
;; 
;; Installation:
;;    add the following lines to your .emacs file:
;;
;;       ; baci
;;       (autoload 'baci-mode "baci" "Enter baci mode." t)
;;       (setq auto-mode-alist (cons '("\\.dat\\'" . baci-mode)
;;                                   auto-mode-alist))
;;       (add-hook 'baci-mode-hook 'turn-on-font-lock)  ;; pbb
;;
;;    copy this file to Emacs site-lisp directory, or individual directory
;;    for instance ~/.site-lisp/baci.
;;
;;       cp baci.el [path to Emacs site-lisp directory] 
;;
;;    If you choose the latter, you need to provide in your .emacs file 
;;    the following line prior to the other lines above
;;
;;       (add-to-list 'load-path "~/.site-lisp/baci")
;;

(defvar baci-mode-hook nil)
(defvar baci-load-hook nil)

(defvar baci-font-lock-keywords
  (list
   '("^//.*" . font-lock-comment-face)  ;; comments
   '("[ ]+//.*" . font-lock-comment-face)  ;; comments
   '("^[-]+[-A-Z0-9 _&]+" . font-lock-reference-face)  ;; sections
   '("[/]+[A-Z0-9 _&]+" . font-lock-type-face)  ;; subsections
;   '("^[A-Z_]+" . font-lock-function-name-face) ;; core keywords
;   '("[A-Z_]+" . font-lock-keyword-face) ;; core keywords
;;   '("[\t]+" . highlight) ;; tabs
;;   '("[ ]+$" . highlight) ;; spaces just before end of line
  )
)

(defun baci-mode ()
  "Major mode for editing baci files."
  (interactive)
  (setq mode-name "baci")
  (setq major-mode 'baci-mode)
  (make-local-variable 'font-lock-defaults)
  (setq font-lock-defaults '(baci-font-lock-keywords nil t))
  (run-hooks 'baci-mode-hook)
)

(provide 'baci-mode)
(run-hooks 'baci-load-hook)


